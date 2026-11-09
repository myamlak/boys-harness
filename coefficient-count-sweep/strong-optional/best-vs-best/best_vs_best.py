#!/usr/bin/env python3
"""Balanced best-vs-best experiment: minimax polynomial against BRASIL rational.

The paper's stored-coefficient comparison is between two pipelines: the
kernel's Chebyshev generator and a rational solver.  That is not a
best-vs-best comparison.  This script supplies the missing polynomial arm: at
each width and target of the corrected BRASIL record it computes the
best-uniform (minimax) polynomial, stores its monomial coefficients in
binary64, evaluates them by Horner on the same grids, and compares the stored
count with the corrected record's rational count.

Matched protocol, cell by cell:

  interval      [x0, x0+W), t = 2(x-x0)/W - 1 in [-1, 1]
  target        rel:1e-13, abs:5e-14, rel:1e-14, as defined in the record
  reference     F_0 at 70 decimal digits (mpmath)
  storage       monomial coefficients rounded to binary64
  evaluation    Horner, worse of fused (math.fma) and separate routes
  grids         2,001 counting nodes, 4,001 check nodes
  guard         rational: min|den| >= 1e-3 max|den|; polynomial: no
                denominator, guard recorded as not applicable

The polynomial solver is an independent Remez exchange in the Chebyshev
basis, iterated in float64 and polished by a 60-digit linear solve at the
converged alternating extrema; the stored coefficients are the exact
Chebyshev-to-monomial conversion at 60 digits, rounded to binary64.  Its
minimax error is checked against an independent linear-programming solve and
against baryrat's BRASIL polynomial mode (which is used only as a quality
reference: its barycentric-to-monomial conversion is numerically unstable at
these degrees and is not used for the stored arm).

The rational arm is the corrected BRASIL record read as deposited.  The new
stored counts are reported in their own columns; the existing pipeline
comparison is not substituted.

Usage:
  python best_vs_best.py --plan
  python best_vs_best.py --out records/best_vs_best.json
  python best_vs_best.py --cell p5 --target rel:1e-13 --max-degree 40
  python best_vs_best.py --selftest
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
import warnings
from pathlib import Path

import mpmath as mp
import numpy as np
from scipy.special import erf

warnings.filterwarnings("ignore")

HERE = Path(__file__).resolve().parent
PAPER = HERE.parent.parent
CORRECTED = PAPER / "results-corrected" / "corrected-counts.json"
CONTROL_SCRIPT = PAPER / "control-60dig" / "control_60dig_brasil.py"
OUT_DEFAULT = HERE / "records" / "best_vs_best.json"
SELFTEST_OUT = HERE / "records" / "selftest.json"
SCHEMA = "c-strong-optional-best-vs-best-1"
CONVERT_DPS = 60
REFERENCE_DPS = 70
GRID_2001 = 2000
GRID_4001 = 4000
DEN_FLOOR = 1e-3
DMAX_DEFAULT = 60


# ---------------------------------------------------------------------------
# Reference and evaluation


def load_control():
    spec = importlib.util.spec_from_file_location("c60_bvb", CONTROL_SCRIPT)
    mod = importlib.util.module_from_spec(spec)
    sys.modules["c60_bvb"] = mod
    spec.loader.exec_module(mod)
    return mod


def f_np(a, W):
    def f(t):
        x = a + W * (np.asarray(t, dtype=float) + 1.0) / 2.0
        return np.sqrt(np.pi) / (2.0 * np.sqrt(x)) * erf(np.sqrt(x))
    return f


def cheb_T(x, d):
    x = np.asarray(x, dtype=float)
    T = np.empty((len(x), d + 1))
    T[:, 0] = 1.0
    if d >= 1:
        T[:, 1] = x
    for j in range(1, d):
        T[:, j + 1] = 2 * x * T[:, j] - T[:, j - 1]
    return T


def cheb_to_mono(c):
    """Chebyshev coefficients to monomial coefficients, at the current dps."""
    d = len(c) - 1
    out = [mp.mpf(0)] * (d + 1)
    Pprev2 = [mp.mpf(0)] * (d + 1)
    Pprev = [mp.mpf(0)] * (d + 1)
    for j in range(d + 1):
        if j == 0:
            P = [mp.mpf(1)]
        elif j == 1:
            P = [mp.mpf(0), mp.mpf(1)]
        else:
            A = [mp.mpf(0)] * (len(Pprev) + 1)
            for i, v in enumerate(Pprev):
                A[i + 1] += 2 * v
            for i, v in enumerate(Pprev2):
                A[i] -= v
            P = A
        for i, v in enumerate(P):
            out[i] += c[j] * v
        Pprev2, Pprev = Pprev, P
    return out


def horner_routes(q, ts, fs_mp):
    """Binary64 fused and separate Horner, both routes; target-kind errors."""
    qd = [float(c) for c in q]
    wr = wa = 0.0
    at_rel = None
    for fused in (True, False):
        for t, fv in zip(ts, fs_mp):
            t = float(t)
            if fused:
                acc = 0.0
                for c in qd[::-1]:
                    acc = math.fma(acc, t, c)
            else:
                acc = 0.0
                for c in qd[::-1]:
                    acc = acc * t + c
            e = abs(acc - float(fv))
            if e > wa:
                wa = e
            r = e / float(fv)
            if r > wr:
                wr, at_rel = r, t
    return {"rel": wr, "abs": wa, "at_rel_t": at_rel}


def exact_error(q, ts, fs_mp, target_kind):
    """The same fit with 60-digit coefficients and arithmetic (no storage floor)."""
    with mp.workdps(CONVERT_DPS):
        qm = [mp.mpf(mp.nstr(c, CONVERT_DPS)) for c in q]
        wr = mp.mpf(0)
        wa = mp.mpf(0)
        for t, fv in zip(ts, fs_mp):
            tm = mp.mpf(float(t))
            acc = mp.mpf(0)
            for c in reversed(qm):
                acc = acc * tm + c
            e = abs(acc - mp.mpf(fv))
            if e > wa:
                wa = e
            r = e / mp.mpf(fv)
            if r > wr:
                wr = r
    return {"rel": float(wr), "abs": float(wa)}


# ---------------------------------------------------------------------------
# Polynomial minimax solver


def remez_float(f, d, n_dense=20001, maxiter=80):
    """Remez exchange in the Chebyshev basis; float64 iteration."""
    k = np.arange(d + 2)
    x = np.sort(np.cos(np.pi * k / (d + 1)))
    c = None
    it = 0
    for it in range(maxiter):
        T = cheb_T(x, d)
        A = np.hstack([T, ((-1.0) ** np.arange(d + 2))[:, None]])
        try:
            sol = np.linalg.solve(A, f(x))
        except np.linalg.LinAlgError:
            return None, x, it, "linear solve failed"
        c = sol[:d + 1]
        xg = np.linspace(-1, 1, n_dense)
        e = f(xg) - cheb_T(xg, d) @ c
        ae = np.abs(e)
        idx = list(np.where((ae[1:-1] >= ae[:-2])
                            & (ae[1:-1] >= ae[2:]))[0] + 1)
        if ae[0] >= ae[1]:
            idx = [0] + idx
        if ae[-1] >= ae[-2]:
            idx = idx + [n_dense - 1]
        comp = []
        for i in sorted(set(idx)):
            if comp and i - comp[-1][-1] <= 2:
                comp[-1].append(i)
            else:
                comp.append([i])
        cand = [max(g, key=lambda i: ae[i]) for g in comp]
        if len(cand) < d + 2:
            return c, x, it, f"only {len(cand)} extrema"
        mags = [ae[i] for i in cand]
        best_j, best_key = 0, -1.0
        for j in range(0, len(cand) - (d + 2) + 1):
            key = min(mags[j:j + d + 2])
            if key > best_key:
                best_j, best_key = j, key
        top = cand[best_j:best_j + d + 2]
        s = np.sign(e[top])
        if not np.all(s[1:] != s[:-1]):
            return c, x, it, "alternation failed"
        xnew = []
        for i in top:
            if i == 0 or i == n_dense - 1:
                xnew.append(xg[i])
            else:
                y0, y1, y2 = e[i - 1], e[i], e[i + 1]
                den = y0 - 2 * y1 + y2
                shift = 0.5 * (y0 - y2) / den if den != 0 else 0.0
                shift = max(-0.5, min(0.5, shift))
                xnew.append(xg[i] + shift * (xg[1] - xg[0]))
        xnew = np.array(sorted(xnew))
        move = float(np.max(np.abs(xnew - x)))
        x = xnew
        if move < 1e-12:
            break
    T = cheb_T(x, d)
    A = np.hstack([T, ((-1.0) ** np.arange(d + 2))[:, None]])
    sol = np.linalg.solve(A, f(x))
    return sol[:d + 1], x, it + 1, None


def f0_mp_of(a, W):
    def f(t):
        x = a + W * (t + 1) / 2
        return mp.sqrt(mp.pi) / (2 * mp.sqrt(x)) * mp.erf(mp.sqrt(x))
    return f


def f0_prime_mp_of(a, W):
    """d/dt F_0(x(t)) with x = a + W(t+1)/2."""
    def fp(t):
        x = a + W * (t + 1) / 2
        dx = (-mp.sqrt(mp.pi) * mp.erf(mp.sqrt(x))
              / (4 * x ** mp.mpf(1.5)) + mp.e ** (-x) / (2 * x))
        return dx * W / 2
    return fp


def pval_cheb(c, t):
    d = len(c) - 1
    b1 = mp.mpf(0)
    b2 = mp.mpf(0)
    for k in range(d, 0, -1):
        b0 = 2 * t * b1 - b2 + c[k]
        b2 = b1
        b1 = b0
    return t * b1 - b2 + c[0]


def pderiv_cheb(c, t):
    d = len(c) - 1
    s = mp.mpf(0)
    if d >= 1:
        u0, u1 = mp.mpf(1), 2 * t
        s += c[1] * u0
        for j in range(2, d + 1):
            u2 = 2 * t * u1 - u0
            s += c[j] * j * u1
            u0, u1 = u1, u2
    return s


def remez_mp(a, W, d, x_ref, f_mp=None, f_prime_mp=None, iters=8,
             ngrid=1001):
    """High-precision Remez exchange at 60 digits.

    Refines the float reference points by solving e'(t)=0 at each extremum,
    re-solves the linear system, and returns the coefficients, the
    equioscillation error E, the amplitudes, and the dense-grid maximum.
    """
    if f_mp is None:
        f_mp = f0_mp_of(a, W)
    if f_prime_mp is None:
        f_prime_mp = f0_prime_mp_of(a, W)
    with mp.workdps(CONVERT_DPS):
        xs = [mp.mpf(float(v)) for v in x_ref]
        c = None
        E = None
        tg = [mp.mpf(-1) + 2 * mp.mpf(k) / (ngrid - 1) for k in range(ngrid)]

        def solve_at(points):
            n = len(points)
            A = mp.matrix(n, d + 2)
            rhs = mp.matrix(n, 1)
            for i, x in enumerate(points):
                T = [mp.mpf(0)] * (d + 1)
                T[0] = mp.mpf(1)
                if d >= 1:
                    T[1] = x
                for j in range(1, d):
                    T[j + 1] = 2 * x * T[j] - T[j - 1]
                for j in range(d + 1):
                    A[i, j] = T[j]
                A[i, d + 1] = (-1) ** i
                rhs[i] = f_mp(x)
            sol = mp.lu_solve(A, rhs)
            return [sol[j] for j in range(d + 1)], sol[d + 1]

        c, E = solve_at(xs)
        best = None                 # (dense, c, E, amps)
        last_dense = None
        for it in range(iters):
            ev = [f_mp(t) - pval_cheb(c, t) for t in tg]
            ae = [abs(v) for v in ev]
            idx = []
            for k in range(1, ngrid - 1):
                if (ae[k] >= ae[k - 1] and ae[k] >= ae[k + 1]
                        and (ae[k] > ae[k - 1] or ae[k] > ae[k + 1])):
                    idx.append(k)
            if ae[0] >= ae[1]:
                idx = [0] + idx
            if ae[-1] >= ae[-2]:
                idx = idx + [ngrid - 1]
            comp = []
            for k in sorted(set(idx)):
                if comp and k - comp[-1][-1] <= 2:
                    comp[-1].append(k)
                else:
                    comp.append([k])
            cand = [max(g, key=lambda k: ae[k]) for g in comp]
            if len(cand) < d + 2:
                break
            mags = [ae[k] for k in cand]
            best_j, best_key = 0, mp.mpf(-1)
            for j in range(len(cand) - (d + 2) + 1):
                key = min(mags[j:j + d + 2])
                if key > best_key:
                    best_j, best_key = j, key
            top = cand[best_j:best_j + d + 2]
            s = [1 if ev[k] > 0 else -1 for k in top]
            if any(s[i] == s[i + 1] for i in range(len(s) - 1)):
                break

            def ederiv(t):
                return f_prime_mp(t) - pderiv_cheb(c, t)

            new = []
            spacing = 2 * mp.mpf(1) / (ngrid - 1)
            for k in top:
                if k == 0:
                    new.append(mp.mpf(-1))
                elif k == ngrid - 1:
                    new.append(mp.mpf(1))
                else:
                    t = tg[k]
                    r = None
                    try:
                        r = mp.findroot(ederiv, t, tol=mp.mpf(10) ** -45)
                    except Exception:
                        r = None
                    try:
                        if r is None or mp.im(r) != 0 or abs(r - t) > 3 * spacing:
                            r = mp.findroot(ederiv, (tg[k - 1], tg[k + 1]),
                                            solver="bisect")
                    except Exception:
                        r = t
                    if r is None or mp.im(r) != 0 or abs(r - t) > 3 * spacing:
                        r = t
                    new.append(mp.re(r))
            move = max(abs(new[i] - xs[i]) for i in range(len(xs)))
            xs = new
            c, E = solve_at(xs)
            last_dense = max(abs(f_mp(t) - pval_cheb(c, t)) for t in tg)
            amps = [f_mp(x) - pval_cheb(c, x) for x in xs]
            if best is None or last_dense < best[0]:
                best = (last_dense, c, E, amps)
            if move < mp.mpf(10) ** -40:
                break
        if best is None:
            last_dense = max(abs(f_mp(t) - pval_cheb(c, t)) for t in tg)
            amps = [f_mp(x) - pval_cheb(c, x) for x in xs]
            best = (last_dense, c, E, amps)
        dense, c, E, amps = best
        valid = False
        if E != 0:
            ratio = dense / abs(E)
            valid = (mp.isfinite(ratio) and mp.mpf("0.5") <= ratio <= 2)
        return c, E, amps, dense, it + 1, bool(valid)


def solve_mp(a, W, d, x_ref, f_mp=None):
    """60-digit linear solve at the converged alternating extrema."""
    if f_mp is None:
        f_mp = f0_mp_of(a, W)
    with mp.workdps(CONVERT_DPS):
        def f(t):
            return f_mp(t)

        n = len(x_ref)
        A = mp.matrix(n, d + 2)
        rhs = mp.matrix(n, 1)
        for i, xv in enumerate(x_ref):
            x = mp.mpf(float(xv))
            T = [mp.mpf(0)] * (d + 1)
            T[0] = mp.mpf(1)
            if d >= 1:
                T[1] = x
            for j in range(1, d):
                T[j + 1] = 2 * x * T[j] - T[j - 1]
            for j in range(d + 1):
                A[i, j] = T[j]
            A[i, d + 1] = (-1) ** i
            rhs[i] = f(x)
        sol = mp.lu_solve(A, rhs)
        c = [sol[j] for j in range(d + 1)]
        E = sol[d + 1]

        def pval(x):
            b1 = mp.mpf(0)
            b2 = mp.mpf(0)
            for k in range(d, 0, -1):
                b0 = 2 * x * b1 - b2 + c[k]
                b2 = b1
                b1 = b0
            return x * b1 - b2 + c[0]

        amps = [f(mp.mpf(float(xv))) - pval(mp.mpf(float(xv)))
                for xv in x_ref]
        return c, E, amps


def equioscillation(amps):
    aa = [abs(v) for v in amps]
    if not aa:
        return None
    mx = max(aa)
    return {
        "n_extrema": len(aa),
        "spread": float((mx - min(aa)) / mx) if mx > 0 else 0.0,
        "amplitudes": [float(v) for v in aa],
    }


# ---------------------------------------------------------------------------
# One cell


def polynomial_arm(cell, control, dmax, verbose=False):
    W = float(cell["width"])
    a = float(cell["_x0"])
    target_kind, target_value = cell["_target_kind"], float(cell["_target_value"])
    mp.mp.dps = REFERENCE_DPS
    ts, xs, fs = control.build_grid(a, W, n=GRID_2001)
    ts4, xs4, fs4 = control.build_grid(a, W, n=GRID_4001)
    f = f_np(a, W)
    scan = []
    accepted = None
    for d in range(0, dmax + 1):
        started = time.perf_counter()
        try:
            c_float, x_ref, iters, err = remez_float(f, d)
        except Exception as ex:  # noqa: BLE001
            scan.append({"degree": d, "stored": d + 1,
                         "exception": f"{type(ex).__name__}: {ex}"})
            continue
        if c_float is None:
            scan.append({"degree": d, "stored": d + 1,
                         "remez_failure": err})
            continue
        with mp.workdps(CONVERT_DPS):
            q_float = cheb_to_mono(c_float)
        ev = horner_routes(q_float, ts, fs)
        err_kind = ev[target_kind]
        meets = err_kind <= target_value
        row = {
            "degree": d, "stored": d + 1,
            "remez_iterations": iters, "remez_note": err,
            "delivered_float_rel": ev["rel"], "delivered_float_abs": ev["abs"],
            "meets_float": bool(meets),
            "seconds": round(time.perf_counter() - started, 4),
        }
        scan.append(row)
        if not meets:
            continue
        # Polish at the first meeting degree; accept only if the polished
        # binary64 storage also meets.
        try:
            c_mp, E, amps = solve_mp(a, W, d, x_ref)
            with mp.workdps(CONVERT_DPS):
                q_mp = cheb_to_mono(c_mp)
            ev2 = horner_routes(q_mp, ts, fs)
            meets2 = ev2[target_kind] <= target_value
            row["polished"] = True
            row["polished_delivered_rel"] = ev2["rel"]
            row["polished_delivered_abs"] = ev2["abs"]
            row["polished_meets"] = bool(meets2)
            row["equioscillation"] = equioscillation(amps)
            row["minimax_abs_error_E"] = float(abs(E))
            if not meets2:
                continue
            accepted = {
                "degree": d,
                "stored": d + 1,
                "coefficients_binary64": [repr(float(c)) for c in q_mp],
                "polished": True,
                "remez_iterations": iters,
                "remez_note": err,
                "minimax_abs_error_E": float(abs(E)),
                "equioscillation": equioscillation(amps),
                "delivered_2001": ev2[target_kind],
                "delivered_rel_2001": ev2["rel"],
                "delivered_abs_2001": ev2["abs"],
                "meets_2001": True,
                "exact_coefficients": exact_error(
                    q_mp, ts, fs, target_kind),
            }
            ev4 = horner_routes(q_mp, ts4, fs4)
            accepted["delivered_4001"] = ev4[target_kind]
            accepted["delivered_rel_4001"] = ev4["rel"]
            accepted["delivered_abs_4001"] = ev4["abs"]
            accepted["meets_4001"] = ev4[target_kind] <= target_value
            accepted["guard"] = None
            accepted["guard_applicable"] = False
            break
        except Exception as ex:  # noqa: BLE001
            row["polish_exception"] = f"{type(ex).__name__}: {ex}"
            continue
    return {"accepted": accepted, "scan": scan}


def rational_arm(cell):
    kind = cell["_target_kind"]
    delivered = cell.get("delivered_2001")
    return {
        "stored": cell["corrected_count"],
        "pair": cell["pair"],
        "route": cell.get("route"),
        "delivered_2001": delivered,
        "delivered_4001": cell.get("delivered_4001"),
        "guard": cell.get("guard"),
        "meets_2001": (
            delivered is not None
            and delivered <= float(cell["_target_value"])),
        "source": "results-corrected/corrected-counts.json, BRASIL 2.1.2",
    }


def run_cells(cells, control, dmax, only_width=None, only_target=None,
              verbose=True):
    out = []
    total = 0
    for cell in cells:
        if only_width and cell["width_name"] != only_width:
            continue
        if only_target and cell["target"] != only_target:
            continue
        total += 1
        t0 = time.perf_counter()
        try:
            poly = polynomial_arm(cell, control, dmax)
        except Exception as ex:  # noqa: BLE001
            poly = {"accepted": None, "scan": [],
                    "exception": f"{type(ex).__name__}: {ex}",
                    "traceback": traceback.format_exc()}
        rat = rational_arm(cell)
        acc = poly.get("accepted")
        row = {
            "width_name": cell["width_name"],
            "width": float(cell["width"]),
            "target": cell["target"],
            "target_kind": cell["_target_kind"],
            "target_value": float(cell["_target_value"]),
            "a": cell["_x0"],
            "rational": rat,
            "polynomial": poly,
            "comparison": None,
        }
        if acc:
            diff = rat["stored"] - acc["stored"]
            row["comparison"] = {
                "rational_stored": rat["stored"],
                "polynomial_stored": acc["stored"],
                "difference_rational_minus_polynomial": diff,
                "case": ("rational_fewer" if diff < 0
                         else "equal" if diff == 0
                         else "polynomial_fewer"),
            }
        out.append(row)
        if verbose:
            if acc:
                print(f"[{total:3d}] {row['width_name']:6s} {row['target']:12s} "
                      f"rat={rat['stored']:2d} poly={acc['stored']:2d} "
                      f"diff={row['comparison']['difference_rational_minus_polynomial']:+d} "
                      f"poly_err={acc['delivered_2001']:.2e} "
                      f"eq={acc['equioscillation']['spread']:.1e} "
                      f"{time.perf_counter()-t0:.1f}s", flush=True)
            else:
                print(f"[{total:3d}] {row['width_name']:6s} {row['target']:12s} "
                      f"rat={rat['stored']:2d} poly=FAIL "
                      f"{time.perf_counter()-t0:.1f}s", flush=True)
    return out


# ---------------------------------------------------------------------------
# Validation against independent minimax solvers


def lp_minimax(a, W, d, xg, f_override=None):
    from scipy.optimize import linprog
    f = f_np(a, W) if f_override is None else f_override
    by = f(xg)
    T = cheb_T(xg, d)
    N = len(xg)
    A = np.vstack([np.hstack([-T, -np.ones((N, 1))]),
                   np.hstack([T, -np.ones((N, 1))])])
    b = np.concatenate([-by, by])
    cobj = np.zeros(d + 2)
    cobj[-1] = 1.0
    res = linprog(cobj, A_ub=A, b_ub=b,
                  bounds=[(None, None)] * (d + 1) + [(0, None)],
                  method="highs")
    if not res.success:
        return None
    return np.asarray(res.x[:d + 1]), float(res.x[-1])


def validate(cells, control, out, dmax=24):
    """Compare Remez on selected cells against BRASIL's polynomial mode and LP.

    BRASIL's internal minimax error is a quality reference; its monomial
    conversion is not used.  The LP is an independent discrete minimax.
    """
    xs = []
    picks = []
    for cell in cells:
        if cell["width_name"] in ("p5", "m7", "p0") and cell["target"] == "rel:1e-13":
            picks.append(cell)
    for cell in picks:
        W = float(cell["width"])
        a = float(cell["_x0"])
        f = f_np(a, W)
        xg = np.linspace(-1, 1, 8001)
        if cell["width_name"] == "p5":
            degrees = [14, 16, 18, 20, 22, 24]
        elif cell["width_name"] == "m7":
            degrees = [1, 2, 3, 4, 5]
        else:
            degrees = [2, 4, 6, 8, 10, 12]
        for d in degrees:
            if d > dmax:
                continue
            rec = {"width_name": cell["width_name"], "target": cell["target"],
                   "degree": d}
            c_float, x_ref, iters, note = remez_float(f, d)
            if c_float is None:
                rec["remez_failure"] = note
                xs.append(rec)
                continue
            try:
                c_mp, E, amps = solve_mp(a, W, d, x_ref)
                rec["remez_minimax_abs"] = float(abs(E))
                rec["equioscillation_spread"] = equioscillation(amps)["spread"]
            except Exception as ex:  # noqa: BLE001
                rec["remez_polish_exception"] = f"{type(ex).__name__}: {ex}"
                xs.append(rec)
                continue
            lp = lp_minimax(a, W, d, xg)
            if lp is not None:
                rec["lp_minimax_abs"] = lp[1]
                rec["lp_over_remez"] = (lp[1] / float(abs(E))
                                        if E != 0 else None)
            try:
                from baryrat import brasil
                r, info = brasil(
                    f, (-1.0, 1.0), (d, 0), info=True, tol=1e-4,
                    maxiter=1000, npi=-30, init_steps=100, max_step_size=0.1)
                rec["brasil_converged"] = bool(info.converged)
                rec["brasil_deviation"] = float(info.deviation)
                rec["brasil_minimax_abs"] = float(info.error)
                rec["brasil_over_remez"] = (
                    float(info.error) / float(abs(E)) if E != 0 else None)
            except Exception as ex:  # noqa: BLE001
                rec["brasil_exception"] = f"{type(ex).__name__}: {ex}"
            xs.append(rec)
            print(f"  validate {cell['width_name']:6s} d={d:2d} "
                  f"remez={rec.get('remez_minimax_abs')} "
                  f"lp_ratio={rec.get('lp_over_remez')} "
                  f"brasil_ratio={rec.get('brasil_over_remez')} "
                  f"eq={rec.get('equioscillation_spread'):.1e}", flush=True)
    return xs


def validate_brasil_conversion(cells, control):
    """Measure BRASIL's polynomial monomial conversion at p5 d=14 and 16.

    BRASIL's internal minimax error is accurate, but the barycentric-to-
    monomial conversion of its polynomial fit is ill-conditioned.  This block
    records the effect so the stored arm is not confused with it.
    """
    cell = next(c for c in cells
                if c["width_name"] == "p5" and c["target"] == "rel:1e-13")
    W = float(cell["width"])
    a = float(cell["_x0"])
    mp.mp.dps = REFERENCE_DPS
    ts, xs, fs = control.build_grid(a, W, n=GRID_2001)
    f = f_np(a, W)
    from baryrat import brasil
    rows = []
    for d in (14, 16):
        r, info = brasil(f, (-1.0, 1.0), (d, 0), info=True, tol=1e-4,
                         maxiter=1000, npi=-30, init_steps=100,
                         max_step_size=0.1)
        P, Q = control.bary_to_monomial(r)
        A = list(P)
        B = list(Q)
        q = [mp.mpf(0)] * max(1, len(A) - len(B) + 1)
        R = A[:]
        while len(R) >= len(B) and not (len(R) == 1 and R[0] == 0):
            dd = len(R) - len(B)
            c = R[-1] / B[-1]
            q[dd] = c
            for i in range(len(B)):
                R[i + dd] -= c * B[i]
            while len(R) > 1 and R[-1] == 0:
                R.pop()
        ev = horner_routes(q, ts, fs)
        rows.append({
            "width_name": "p5", "width": W, "degree": d,
            "brasil_internal_error": float(info.error),
            "brasil_converged": bool(info.converged),
            "max_abs_monomial_coefficient": float(
                max(abs(float(c)) for c in q)),
            "delivered_rel": ev["rel"],
            "delivered_abs": ev["abs"],
            "note": "BRASIL polynomial mode; barycentric-to-monomial "
                    "conversion at 70 digits then binary64",
        })
    return rows


def selftest():
    cases = []

    def check(name, f, a, W, d, f_mp=None, check_lp=True):
        c, x, it, note = remez_float(f, d)
        rec = {"name": name, "degree": d, "note": note}
        if c is None:
            rec["ok"] = False
            cases.append(rec)
            return
        c_mp, E, amps = solve_mp(a, W, d, x, f_mp=f_mp)
        rec["minimax_abs"] = float(abs(E))
        rec["equioscillation_spread"] = equioscillation(amps)["spread"]
        # The claimed minimax error must be attained on the interval.
        with mp.workdps(CONVERT_DPS):
            q = cheb_to_mono(c_mp)
        qd = [float(v) for v in q]
        xg = np.linspace(-1, 1, 8001)
        pv = np.polyval(qd[::-1], xg)
        dense_err = float(np.max(np.abs(f(xg) - pv)))
        rec["dense_max_error"] = dense_err
        rec["dense_over_E"] = dense_err / float(abs(E)) if E != 0 else None
        lp = lp_minimax(a, W, d, xg, f_override=f)
        if lp is not None:
            rec["lp_minimax_abs"] = lp[1]
            rec["lp_over_remez"] = lp[1] / float(abs(E)) if E != 0 else None
        ok = (rec["equioscillation_spread"] < 1e-8
              and dense_err <= 2.0 * float(abs(E)) + 1e-13)
        if check_lp and lp is not None and float(abs(E)) > 1e-9:
            ok = ok and abs(lp[1] / float(abs(E)) - 1) < 1e-3
        rec["ok"] = ok
        cases.append(rec)

    def fexp(t):
        return np.exp(np.asarray(t, dtype=float))

    check("exp", fexp, 0.0, 2.0, 5, f_mp=lambda t: mp.exp(t))
    check("exp", fexp, 0.0, 2.0, 12, f_mp=lambda t: mp.exp(t), check_lp=False)
    check("F0-m7", f_np(11.899848152108484, 0.0078125),
          11.899848152108484, 0.0078125, 3,
          f_mp=f0_mp_of(11.899848152108484, 0.0078125), check_lp=False)
    check("F0-p5", f_np(11.899848152108484, 32.0),
          11.899848152108484, 32.0, 16,
          f_mp=f0_mp_of(11.899848152108484, 32.0), check_lp=False)
    passed = sum(1 for c in cases if c["ok"])
    return {"cases": cases, "passed": passed, "total": len(cases),
            "all_pass": passed == len(cases)}


# ---------------------------------------------------------------------------


def refine_cell(row, control, max_steps=4, marginal=1.5):
    """High-precision recheck of the degrees below the accepted one.

    The record's accepted degree is kept unless a lower degree meets at high
    precision.  At each lower degree the 60-digit linear polish is evaluated;
    when its delivered error is within `marginal` times the target, the
    high-precision exchange also runs, and the better of the two is used.
    """
    W = float(row["width"])
    a = float(row["a"])
    kind = row["target_kind"]
    tgt = float(row["target_value"])
    d_acc = row["polynomial"]["accepted"]["degree"]
    mp.mp.dps = REFERENCE_DPS
    ts, xs, fs = control.build_grid(a, W, n=GRID_2001)
    ts4, xs4, fs4 = control.build_grid(a, W, n=GRID_4001)
    f = f_np(a, W)
    f_mp = f0_mp_of(a, W)
    checked = []
    refined = None
    d = d_acc - 1
    while d >= 0 and d >= d_acc - max_steps:
        started = time.perf_counter()
        c_float, x_ref, iters, note = remez_float(f, d)
        if c_float is None:
            checked.append({"degree": d, "remez_float_failure": note})
            break
        cands = []

        def evaluate(c_mp, method, valid, dense=None, E=None, amps=None,
                     it2=None):
            with mp.workdps(CONVERT_DPS):
                q = cheb_to_mono(c_mp)
            if any(mp.im(v) != 0 for v in q):
                raise ValueError("complex monomial coefficient")
            ev = horner_routes(q, ts, fs)
            rec = {
                "degree": d,
                "stored": d + 1,
                "method": method,
                "mp_exchange_valid": bool(valid),
                "delivered_2001": ev[kind],
                "delivered_rel_2001": ev["rel"],
                "delivered_abs_2001": ev["abs"],
                "meets": bool(ev[kind] <= tgt),
                "minimax_abs_error_E": (float(abs(E)) if E is not None else None),
                "dense_max_error": (float(dense) if dense is not None else None),
                "dense_over_E": (float(dense / abs(E))
                                 if (dense is not None and E) else None),
                "equioscillation_spread": (equioscillation(amps)["spread"]
                                           if amps is not None else None),
                "remez_mp_iterations": it2,
                "coefficients_binary64": [repr(float(v)) for v in q],
                "coefficients_60": [mp.nstr(v, CONVERT_DPS) for v in q],
                "seconds": round(time.perf_counter() - started, 3),
            }
            rec["_q"] = q
            return rec

        try:
            c1, E1, amps1 = solve_mp(a, W, d, x_ref)
            rec1 = evaluate(c1, "60-digit linear polish", False,
                            E=E1, amps=amps1)
            cands.append(rec1)
        except Exception as ex:  # noqa: BLE001
            checked.append({"degree": d, "linear_polish_exception":
                            f"{type(ex).__name__}: {ex}"})

        delivered1 = min((c["delivered_2001"] for c in cands), default=math.inf)
        if delivered1 <= marginal * tgt:
            try:
                c2, E2, amps2, dense2, it2, valid = remez_mp(
                    a, W, d, x_ref)
                if valid:
                    rec2 = evaluate(c2, "high-precision Remez exchange", True,
                                    dense=dense2, E=E2, amps=amps2, it2=it2)
                    cands.append(rec2)
                else:
                    checked.append({"degree": d,
                                    "mp_exchange_note": "invalid dense/E"})
            except Exception as ex:  # noqa: BLE001
                checked.append({"degree": d, "mp_exchange_exception":
                                f"{type(ex).__name__}: {ex}"})

        if not cands:
            break
        best = min(cands, key=lambda c: c["delivered_2001"])
        q_best = best.pop("_q")
        for c in cands:
            c.pop("_q", None)
        eval4 = horner_routes(q_best, ts4, fs4)
        best["delivered_4001"] = eval4[kind]
        best["meets_4001"] = bool(eval4[kind] <= tgt)
        best["exact_coefficients"] = exact_error(q_best, ts, fs, kind)
        checked.append(best)
        if best["meets"]:
            refined = dict(best)
            d -= 1
        else:
            break
    return {"checked": checked, "refined_accepted": refined}


def recheck(record_path, out_path, control):
    record = json.loads(Path(record_path).read_text())
    rows = record["cells"]
    changed = []
    for n, row in enumerate(rows, 1):
        if not row.get("polynomial", {}).get("accepted"):
            continue
        old_stored = row["polynomial"]["accepted"]["stored"]
        rr = refine_cell(row, control)
        rr["original_polynomial_stored"] = old_stored
        rr["seconds_total"] = round(sum(
            c.get("seconds", 0.0) for c in rr.get("checked", [])), 3)
        row["high_precision_recheck"] = rr
        best = rr.get("refined_accepted")
        if best is not None:
            old = row["polynomial"]["accepted"]
            merged = dict(old)
            merged.update({
                "degree": best["degree"],
                "stored": best["stored"],
                "coefficients_binary64": best["coefficients_binary64"],
                "delivered_2001": best["delivered_2001"],
                "delivered_rel_2001": best["delivered_rel_2001"],
                "delivered_abs_2001": best["delivered_abs_2001"],
                "delivered_4001": best.get("delivered_4001"),
                "meets_2001": True,
                "meets_4001": best.get("meets_4001"),
                "minimax_abs_error_E": best["minimax_abs_error_E"],
                "equioscillation": {
                    "n_extrema": None,
                    "spread": best["equioscillation_spread"],
                },
                "exact_coefficients": best.get("exact_coefficients"),
                "solver": "high-precision Remez exchange (60-digit polish)",
                "dense_max_error": best["dense_max_error"],
                "dense_over_E": best["dense_over_E"],
                "guard": None,
                "guard_applicable": False,
            })
            row["polynomial"]["accepted"] = merged
        rat = row["rational"]
        acc = row["polynomial"]["accepted"]
        diff = rat["stored"] - acc["stored"]
        row["comparison"] = {
            "rational_stored": rat["stored"],
            "polynomial_stored": acc["stored"],
            "difference_rational_minus_polynomial": diff,
            "case": ("rational_fewer" if diff < 0
                     else "equal" if diff == 0
                     else "polynomial_fewer"),
        }
        checked = rr.get("checked", [])
        new_stored = acc["stored"]
        if old_stored != new_stored:
            changed.append({"width_name": row["width_name"],
                            "target": row["target"],
                            "original_polynomial_stored": old_stored,
                            "rechecked_polynomial_stored": new_stored})
        print(f"[{n:3d}] {row['width_name']:6s} {row['target']:12s} "
              f"rat={rat['stored']:2d} poly={acc['stored']:2d} "
              f"diff={diff:+d} rechecked="
              f"{[c.get('degree') for c in checked]} "
              f"meets={[c.get('meets') for c in checked]}", flush=True)

    ok = [r for r in rows if r.get("comparison")]
    cases = {"rational_fewer": 0, "equal": 0, "polynomial_fewer": 0}
    for r in ok:
        cases[r["comparison"]["case"]] += 1
    diffs = [r["comparison"]["difference_rational_minus_polynomial"] for r in ok]
    record["summary"]["cases"] = cases
    record["summary"]["difference_rational_minus_polynomial"] = {
        "min": min(diffs) if diffs else None,
        "max": max(diffs) if diffs else None,
        "mean": (sum(diffs) / len(diffs)) if diffs else None,
        "counts": {str(v): diffs.count(v) for v in sorted(set(diffs))},
    }
    for tgt in ("rel:1e-13", "abs:5e-14", "rel:1e-14"):
        sub = [r for r in ok if r["target"] == tgt]
        record["summary"]["by_target"][tgt] = {
            "cells": len(sub),
            "rational_fewer": sum(1 for r in sub
                                  if r["comparison"]["case"] == "rational_fewer"),
            "equal": sum(1 for r in sub
                         if r["comparison"]["case"] == "equal"),
            "polynomial_fewer": sum(1 for r in sub
                                    if r["comparison"]["case"] == "polynomial_fewer"),
            "mean_difference": (
                sum(r["comparison"]["difference_rational_minus_polynomial"]
                    for r in sub) / len(sub)) if sub else None,
        }
    for w in [c["width_name"] for c in rows[:15]]:
        sub = [r for r in rows if r["width_name"] == w]
        if sub:
            record["summary"]["by_width"][w] = [
                {"target": r["target"],
                 "rational": r["rational"]["stored"],
                 "polynomial": r["polynomial"]["accepted"]["stored"],
                 "difference": r["comparison"][
                     "difference_rational_minus_polynomial"]}
                for r in sub]
    record["summary"]["recheck_seconds_total"] = round(sum(
        r.get("high_precision_recheck", {}).get("seconds_total", 0.0)
        for r in rows), 2)
    record["summary"]["counts_changed_by_recheck"] = changed
    Path(out_path).write_text(json.dumps(record, indent=2) + "\n")
    print(json.dumps(record["summary"], indent=2))
    print(f"record {out_path}")
    return record


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def load_cells():
    doc = json.loads(CORRECTED.read_text())
    x0 = float(doc["problem"]["x0"])
    widths = {k: float(v) for k, v in doc["problem"]["widths"].items()}
    cells = []
    for c in doc["cells"]:
        kind, val = c["target"].split(":")
        cell = dict(c)
        cell["_x0"] = x0
        cell["_width_value"] = widths[c["width_name"]]
        cell["_target_kind"] = "rel" if kind == "rel" else "abs"
        cell["_target_value"] = float(val)
        cells.append(cell)
    return doc, cells


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", type=Path, default=OUT_DEFAULT)
    ap.add_argument("--plan", action="store_true")
    ap.add_argument("--selftest", action="store_true")
    ap.add_argument("--cell", default=None)
    ap.add_argument("--target", default=None)
    ap.add_argument("--max-degree", type=int, default=DMAX_DEFAULT)
    ap.add_argument("--no-validate", action="store_true")
    ap.add_argument("--recheck", action="store_true",
                    help="high-precision recheck of the accepted and below "
                         "degrees on an existing record")
    ap.add_argument("--recheck-input", type=Path, default=OUT_DEFAULT)
    args = ap.parse_args(argv)

    if args.selftest:
        st = selftest()
        SELFTEST_OUT.parent.mkdir(parents=True, exist_ok=True)
        SELFTEST_OUT.write_text(json.dumps(st, indent=2) + "\n")
        for c in st["cases"]:
            print(f"  {'ok  ' if c['ok'] else 'FAIL'} {c['name']} d={c['degree']} "
                  f"E={c.get('minimax_abs')} eq={c.get('equioscillation_spread')} "
                  f"lp_ratio={c.get('lp_over_remez')} note={c.get('note')}")
        print(f"selftest: {st['passed']}/{st['total']} passed; "
              f"record {SELFTEST_OUT}")
        return 0 if st["all_pass"] else 1

    doc, cells = load_cells()
    if args.recheck:
        control = load_control()
        recheck(args.recheck_input, args.out, control)
        return 0
    if args.plan:
        for c in cells:
            print(f"{c['width_name']:6s} {c['width']:>10} {c['target']:12s} "
                  f"rational={c['corrected_count']:2d} pair={c['pair']}")
        print(f"cells: {len(cells)}")
        return 0

    control = load_control()
    started = time.perf_counter()
    rows = run_cells(cells, control, args.max_degree,
                     only_width=args.cell, only_target=args.target)
    validation = []
    if not args.no_validate:
        print("validation:", flush=True)
        validation = validate(cells, control, args.out)
        conversion = validate_brasil_conversion(cells, control)
    else:
        conversion = []

    ok = [r for r in rows if r["comparison"]]
    cases = {"rational_fewer": 0, "equal": 0, "polynomial_fewer": 0}
    for r in ok:
        cases[r["comparison"]["case"]] += 1
    diffs = [r["comparison"]["difference_rational_minus_polynomial"] for r in ok]
    summary = {
        "cells": len(rows),
        "cells_with_accepted_polynomial": len(ok),
        "cells_polynomial_failed": len(rows) - len(ok),
        "cases": cases,
        "difference_rational_minus_polynomial": {
            "min": min(diffs) if diffs else None,
            "max": max(diffs) if diffs else None,
            "mean": (sum(diffs) / len(diffs)) if diffs else None,
            "counts": {str(v): diffs.count(v) for v in sorted(set(diffs))},
        },
        "by_target": {},
        "by_width": {},
        "seconds": round(time.perf_counter() - started, 2),
    }
    for tgt in ("rel:1e-13", "abs:5e-14", "rel:1e-14"):
        sub = [r for r in ok if r["target"] == tgt]
        summary["by_target"][tgt] = {
            "cells": len(sub),
            "rational_fewer": sum(1 for r in sub
                                  if r["comparison"]["case"] == "rational_fewer"),
            "equal": sum(1 for r in sub
                         if r["comparison"]["case"] == "equal"),
            "polynomial_fewer": sum(1 for r in sub
                                    if r["comparison"]["case"] == "polynomial_fewer"),
            "mean_difference": (
                sum(r["comparison"]["difference_rational_minus_polynomial"]
                    for r in sub) / len(sub)) if sub else None,
        }
    for w in [c["width_name"] for c in cells[:15]]:
        sub = [r for r in ok if r["width_name"] == w]
        if sub:
            summary["by_width"][w] = [
                {"target": r["target"],
                 "rational": r["rational"]["stored"],
                 "polynomial": r["polynomial"]["accepted"]["stored"],
                 "difference": r["comparison"][
                     "difference_rational_minus_polynomial"]}
                for r in sub]

    record = {
        "schema": SCHEMA,
        "generated_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "versions": {
            "python": sys.version.split()[0],
            "mpmath": mp.__version__,
            "numpy": np.__version__,
            "scipy": __import__("scipy").__version__,
            "baryrat": __import__("baryrat").__version__,
        },
        "script": {"path": str(Path(__file__).resolve()),
                   "sha256": sha256(Path(__file__).resolve())},
        "inputs": {"corrected_record": {
            "path": str(CORRECTED), "sha256": sha256(CORRECTED)}},
        "protocol": {
            "interval": "t = 2(x-x0)/W - 1 in [-1,1]",
            "reference": "F_0 at 70 decimal digits (mpmath)",
            "storage": "monomial coefficients rounded to binary64",
            "evaluation": "Horner, worse of fused (math.fma) and separate "
                          "routes, denominator guard min|den| >= 1e-3 max|den|",
            "counting_grid": GRID_2001 + 1,
            "check_grid": GRID_4001 + 1,
            "polynomial_solver": "Remez exchange in the Chebyshev basis "
                                 "(float64) + 60-digit linear polish; exact "
                                 "Chebyshev-to-monomial conversion at 60 digits",
            "rational_arm": "the corrected BRASIL record as deposited",
            "guard": "rational only; the polynomial has no denominator",
            "search": "degrees from 0 upward; the stored count is the first "
                      "degree whose binary64-delivered error meets the target",
        },
        "validation": validation,
        "validation_brasil_polynomial_conversion": conversion,
        "cells": rows,
        "summary": summary,
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(record, indent=2) + "\n")
    print(json.dumps(summary, indent=2))
    print(f"record {args.out}")
    return 0 if summary["cells_polynomial_failed"] == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
