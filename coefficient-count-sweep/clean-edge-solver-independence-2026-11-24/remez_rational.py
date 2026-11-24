"""Independent high-precision rational minimax solver for the clean-edge cells.

Algorithm family: classical rational Remez exchange (differential correction
with a nonlinear Newton polish), implemented from the algorithm at 70 decimal
digits with mpmath.  This is independent of BRASIL/baryrat (barycentric Lawson
iteration) and of the record's binary64 kernel Remez/Lawson solver.

Problem solved: for t in [-1, 1] and a positive weight w(t), find P_m/Q_k,
Q_k(0)=1, minimising

    max_t  w(t) |P(t)/Q(t) - f(t)|.

The clean-edge control uses w = 1/f, so the minimised quantity is the relative
error, exactly the control's acceptance metric.  The acceptance and the
reported e60 are evaluated on the control's 2,001-node uniform-t grid.

Reference set: N = m + k + 2 points.  Alternation theorem: if the weighted
error takes values of equal magnitude |E| with alternating signs at N reference
points, and the maximum over the interval is |E|, the rational is the best
approximation.  The solver reports a fit as converged when the reference
errors alternate with equal magnitude and the grid maximum is within CONV_TOL
of that magnitude; the certificate, not the iteration history, is the
convergence criterion.

Two phases:
  * differential correction -- the classical linearised Remez step; it solves
        w_i (P(t_i) - f_i Q(t_i)) - s_i E = 0,   s_i = (-1)^i,
    for the coefficients and the level E (the dropped term is second order);
  * Newton polish -- the exact reference equations
        w_i (P(t_i)/Q(t_i) - f_i) - s_i E = 0,
    solved by Newton at 70 digits, which makes the reference errors exactly
    equioscillating.

Reference points are continuous: each exchange locates the local extrema of
the error by parabolic interpolation on the dense grid, so the reference is
not restricted to grid nodes.  The exchange keeps the N alternating extrema
with the largest minimum magnitude (the de la Vallee Poussin bound).

Known-answer test (`--test`): the minimax polynomial of degree n-1 for t^n on
[-1,1] is t^n - T_n/2^(n-1), with error 2^(1-n).

Usage:
  python -I remez_rational.py --test
  python -I remez_rational.py --a 62 --width 166.72481653881133 --m 6 --k 8
"""

from __future__ import annotations

import argparse
import json
import math
import time
from pathlib import Path

import mpmath as mp

DPS = 70                    # solver and F_0 reference precision
EVAL_DPS = 60               # the control's evaluation precision
K_MAX = 12
DEN_FLOOR = 1e-3
CONV_TOL = mp.mpf("1e-10")  # grid max / reference level - 1 at convergence
DC_ITERS = 30
MAX_OUTER = 100
NEWTON_ITERS = 40
NEWTON_RES_TOL = mp.mpf("1e-38")
NEWTON_STEP_TOL = mp.mpf("1e-45")

mp.mp.dps = DPS


# ---------------------------------------------------------------------------
# The clean-edge function and grid


def f0_mp(x):
    """F_0(x) = sqrt(pi)/(2 sqrt(x)) erf(sqrt(x)) at the current precision."""
    return mp.sqrt(mp.pi) / (2 * mp.sqrt(x)) * mp.erf(mp.sqrt(x))


def build_grid(a, W, n=2000):
    """The control's grid: binary64 t and x, F_0 at 70 decimal digits."""
    ts = []
    tsm = []
    fs = []
    for i in range(n + 1):
        t = 2.0 * i / n - 1.0               # binary64, as in the control
        x = a + W * (t + 1.0) / 2.0         # binary64, as in the control
        ts.append(t)
        tsm.append(mp.mpf(t))
        fs.append(f0_mp(mp.mpf(float(x))))
    return ts, tsm, fs


class Problem:
    """f(t) and the weight w(t) at arbitrary t."""

    def __init__(self, eval_f, weight):
        self.eval_f = eval_f
        self.weight = weight

    def f(self, t):
        return self.eval_f(t)

    def w(self, t):
        return self.weight(t, self.eval_f(t))


def problem_f0(a, W):
    def eval_f(t):
        x = mp.mpf(a) + mp.mpf(W) * (t + 1) / 2
        return f0_mp(x)
    return Problem(eval_f, lambda t, f: 1 / f)


def problem_power(n_deg):
    return Problem(lambda t: t ** n_deg, lambda t, f: mp.mpf(1))


# ---------------------------------------------------------------------------
# Polynomial and error evaluation


def poly_eval(cs, t):
    acc = mp.mpf(0)
    for c in reversed(cs):
        acc = acc * t + c
    return acc


def profile(P, Q, tsm, fs, ws):
    """Signed weighted error and Q values at every grid node."""
    errs = []
    qvals = []
    for t, f, w in zip(tsm, fs, ws):
        q = poly_eval(Q, t)
        p = poly_eval(P, t)
        errs.append(w * (p / q - f))
        qvals.append(q)
    return errs, qvals


def error_at(P, Q, t, prob):
    q = poly_eval(Q, t)
    p = poly_eval(P, t)
    f = prob.f(t)
    w = prob.weight(t, f)
    return w * (p / q - f), q


def guard_ratio(qvals):
    dmin = min(abs(q) for q in qvals)
    dmax = max(abs(q) for q in qvals)
    if dmax == 0:
        return mp.mpf(0)
    return dmin / dmax


# ---------------------------------------------------------------------------
# Reference initialisation and exchange


def chebyshev_reference(N):
    """The N extrema of T_{N-1} on [-1,1], as exact mpf values."""
    return [-mp.cos(mp.pi * i / (N - 1)) for i in range(N)]


def parabola_vertex(tm, y0, y1, y2):
    """Vertex of the parabola through (-h, y0), (0, y1), (h, y2)."""
    h = tm[1] - tm[0]
    den = y0 - 2 * y1 + y2
    if den == 0:
        return tm[1]
    off = h * (y0 - y2) / (2 * den)
    if abs(off) > 2 * h:
        return tm[1]
    return tm[1] + off


def refine_extrema(P, Q, ts, errs, prob):
    """Refined local extrema of the signed weighted error, alternating.

    Returns a list of (t, e) pairs merged to strict sign alternation.
    """
    n = len(errs) - 1
    cand = [0]
    for i in range(1, n):
        a, b, c = abs(errs[i - 1]), abs(errs[i]), abs(errs[i + 1])
        if b >= a and b >= c and (b > a or b > c):
            cand.append(i)
    cand.append(n)
    refined = []
    for i in cand:
        if i == 0:
            t = mp.mpf(ts[0])
        elif i == n:
            t = mp.mpf(ts[n])
        else:
            s = 1 if errs[i] > 0 else -1
            # signed error, maximised at the extremum
            y0, y1, y2 = s * errs[i - 1], s * errs[i], s * errs[i + 1]
            t = parabola_vertex(ts[i - 1:i + 2], y0, y1, y2)
            t = mp.mpf(max(ts[i - 1], min(ts[i + 1], t)))
        e, _ = error_at(P, Q, t, prob)
        refined.append((t, e))
    # merge same-sign neighbours, keeping the larger magnitude
    merged = []
    for t, e in refined:
        if e == 0:
            continue
        s = 1 if e > 0 else -1
        if merged and (1 if merged[-1][1] > 0 else -1) == s:
            if abs(e) > abs(merged[-1][1]):
                merged[-1] = (t, e)
        else:
            merged.append((t, e))
    return merged


def choose_window(merged, N):
    """The N alternating extrema with the largest minimum magnitude.

    The merged extrema alternate in sign, so an alternating subsequence is any
    selection whose consecutive indices differ by an odd number.  Dynamic
    programming finds the subsequence of length N maximising the smallest
    magnitude -- the strongest de la Vallee Poussin lower bound.
    """
    M = len(merged)
    if M < N:
        return None
    mags = [abs(e) for _, e in merged]
    # dp[c][j]: best minimum magnitude for c+1 points ending at j
    dp = [[None] * M for _ in range(N)]
    prev = [[None] * M for _ in range(N)]
    for j in range(M):
        dp[0][j] = mags[j]
    for c in range(1, N):
        for j in range(M):
            best = None
            bi = None
            for i in range(j):
                if (j - i) % 2 == 1 and dp[c - 1][i] is not None:
                    val = min(dp[c - 1][i], mags[j])
                    if best is None or val > best:
                        best, bi = val, i
            dp[c][j] = best
            prev[c][j] = bi
    candidates = [(dp[N - 1][j], j) for j in range(M)
                  if dp[N - 1][j] is not None]
    if not candidates:
        return None
    _, j = max(candidates)
    seq = [j]
    for c in range(N - 1, 0, -1):
        j = prev[c][j]
        seq.append(j)
    seq.reverse()
    return [merged[i][0] for i in seq]


def maehly_exchange(merged, ref):
    """Remez second algorithm: move each reference point to the local
    extremum of the error in its Voronoi cell (between the midpoints of the
    neighbouring reference points).  Slower than a global exchange but more
    conservative, which keeps the iteration out of spurious near-pole basins.
    """
    N = len(ref)
    mids = [(ref[j] + ref[j + 1]) / 2 for j in range(N - 1)]
    new = []
    for j in range(N):
        lo = mids[j - 1] if j > 0 else None
        hi = mids[j] if j < N - 1 else None
        cands = [(t, e) for t, e in merged
                 if (lo is None or t >= lo) and (hi is None or t <= hi)]
        if not cands:
            new.append(ref[j])
            continue
        new.append(min(cands, key=lambda x: abs(x[0] - ref[j]))[0])
    # keep the reference strictly increasing; an inversion is a failed move
    for j in range(1, N):
        if new[j] <= new[j - 1]:
            new[j] = ref[j]
    return new


def alternating_signs_ok(values):
    signs = []
    for e in values:
        if e == 0:
            return False
        signs.append(1 if e > 0 else -1)
    return all(signs[i] == -signs[i - 1] for i in range(1, len(signs)))


# ---------------------------------------------------------------------------
# Reference solves


def dc_solve(m, k, ref, prob):
    """One differential-correction step on the reference set.

    Rows:  w_i t_i^j  (p_j),  -w_i f_i t_i^j  (q_j),  -s_i  (E)  =  w_i f_i.
    """
    A = []
    b = []
    for i, t in enumerate(ref):
        f = prob.f(t)
        w = prob.weight(t, f)
        s = 1 if i % 2 == 0 else -1
        row = []
        tp = mp.mpf(1)
        for _ in range(m + 1):
            row.append(w * tp)
            tp *= t
        tp = mp.mpf(1)
        for _ in range(k):
            tp *= t
            row.append(-w * f * tp)
        row.append(mp.mpf(-s))
        A.append(row)
        b.append(w * f)
    u = mp.lu_solve(mp.matrix(A), mp.matrix(b))
    P = [u[j] for j in range(m + 1)]
    Q = [mp.mpf(1)] + [u[m + 1 + j] for j in range(k)]
    E = u[m + 1 + k]
    return P, Q, E


def newton_solve(m, k, ref, prob, P, Q, E):
    """Newton on the exact weighted equioscillation equations.

    R_i = w_i (P(t_i)/Q(t_i) - f_i) - s_i E = 0.
    """
    P = list(P)
    Q = list(Q)
    if len(Q) < k + 1:
        Q += [mp.mpf(0)] * (k + 1 - len(Q))
    E = mp.mpf(E)
    last_res = None
    for _ in range(NEWTON_ITERS):
        J = []
        R = []
        for i, t in enumerate(ref):
            f = prob.f(t)
            w = prob.weight(t, f)
            s = 1 if i % 2 == 0 else -1
            Pi = poly_eval(P, t)
            Qi = poly_eval(Q, t)
            if Qi == 0:
                return P, Q, E, False, "Q(t_i) = 0", last_res
            R.append(w * (Pi / Qi - f) - s * E)
            row = []
            tp = mp.mpf(1)
            for _ in range(m + 1):
                row.append(w * tp)
                tp *= t
            tp = mp.mpf(1)
            for _ in range(k):
                tp *= t
                row.append(-w * f * tp - s * E * tp)
            row.append(mp.mpf(-s) * Qi)
            J.append(row)
        last_res = max(abs(r) for r in R)
        if last_res <= NEWTON_RES_TOL:
            return P, Q, E, True, f"residual {mp.nstr(last_res, 5)}", last_res
        try:
            delta = mp.lu_solve(mp.matrix(J), mp.matrix([-r for r in R]))
        except (ZeroDivisionError, ValueError) as ex:
            return P, Q, E, False, f"singular Jacobian: {ex!r}", last_res
        step = max(abs(delta[j]) for j in range(len(delta)))
        if step <= NEWTON_STEP_TOL:
            return P, Q, E, True, f"step {mp.nstr(step, 5)}", last_res

        # backtracking line search on the exact residual norm
        alpha = mp.mpf(1)
        accepted = False
        for _ in range(14):
            Pt = [P[j] + alpha * delta[j] for j in range(m + 1)]
            Qt = [mp.mpf(1)] + [Q[1 + j] + alpha * delta[m + 1 + j]
                                for j in range(k)]
            Et = E + alpha * delta[m + 1 + k]
            new_res = _newton_residual(m, k, ref, prob, Pt, Qt, Et)
            if new_res < last_res:
                P, Q, E = Pt, Qt, Et
                accepted = True
                break
            alpha /= 2
        if not accepted:
            return P, Q, E, False, "Newton line search failed", last_res
    return P, Q, E, False, "maximum Newton iterations", last_res


def _newton_residual(m, k, ref, prob, P, Q, E):
    R = []
    for i, t in enumerate(ref):
        f = prob.f(t)
        w = prob.weight(t, f)
        s = 1 if i % 2 == 0 else -1
        R.append(w * (poly_eval(P, t) / poly_eval(Q, t) - f) - s * E)
    return max(abs(r) for r in R)


# ---------------------------------------------------------------------------
# The solver


def solve_pair(m, k, tsm, fs, prob, target,
               dc_iters=DC_ITERS, max_outer=MAX_OUTER, verbose=False,
               ref0=None, init=None, exchange_mode="global"):
    """Minimise the weighted error for one [m/k] on the supplied grid.

    ref0 supplies an initial reference set (default: Chebyshev extrema).
    init=(P,Q,E) warm-starts the fit; with dc_iters=0 the first solve is the
    exact Newton on ref0, which is how degree continuation is run.
    """
    started = time.perf_counter()
    ws = [prob.weight(t, f) for t, f in zip(tsm, fs)]
    N = m + k + 2
    ref = list(ref0) if ref0 is not None else chebyshev_reference(N)
    if init is None:
        P, Q, E = dc_solve(m, k, ref, prob)
    else:
        P, Q, E = list(init[0]), list(init[1]), mp.mpf(init[2])
        if len(Q) < k + 1:
            Q += [mp.mpf(0)] * (k + 1 - len(Q))

    best = None          # best admissible fit by grid max
    cert = None          # equioscillation certificate
    newton_used = False
    newton_ok = False
    newton_note = None
    newton_res = None
    history = []
    converged = False
    stall = 0
    phase = "dc"

    def consider(P, Q, E, ref, it, phase):
        nonlocal best, cert, converged
        errs, qvals = profile(P, Q, tsm, fs, ws)
        g = guard_ratio(qvals)
        gm = max(abs(e) for e in errs)
        ref_errs = [error_at(P, Q, t, prob)[0] for t in ref]
        signs_ok = alternating_signs_ok(ref_errs)
        ref_min = min(abs(e) for e in ref_errs)
        ratio = gm / ref_min if ref_min > 0 else mp.inf
        entry = {"it": it, "phase": phase, "grid_max": gm,
                 "guard": g, "signs_ok": signs_ok, "ratio": ratio,
                 "E": E, "ref_min": ref_min}
        history.append(entry)
        if g >= DEN_FLOOR and (best is None or gm < best[0]["grid_max"]):
            best = (entry, P, Q, E, ref, errs, qvals)
        if signs_ok and ratio <= 1 + CONV_TOL:
            cert = (entry, P, Q, E, ref, errs, qvals)
            converged = True
        return errs, qvals, gm, g, signs_ok, ratio

    errs, qvals, gm, g, signs_ok, ratio = consider(P, Q, E, ref, 0, "dc")
    if verbose:
        print(f"    [start] N={N} grid_max={mp.nstr(gm, 8)} "
              f"guard={mp.nstr(g, 4)} ratio={mp.nstr(ratio, 6)}", flush=True)

    for it in range(1, max_outer + 1):
        if converged:
            break
        merged = refine_extrema(P, Q, tsm, errs, prob)
        if exchange_mode == "maehly":
            new_ref = maehly_exchange(merged, ref)
        else:
            new_ref = choose_window(merged, N)
        if verbose:
            print(f"    it={it:3d} {phase:6s} grid_max={mp.nstr(gm, 8)} "
                  f"ratio={mp.nstr(ratio, 6)} guard={mp.nstr(g, 4)} "
                  f"extrema={len(merged)} signs={signs_ok}", flush=True)
        if new_ref is None:
            # fewer than N alternating extrema: try Newton once, then stop
            if phase == "dc" and not newton_used:
                phase = "newton"
            else:
                break
            new_ref = ref
        if new_ref == ref:
            stall += 1
            if stall > 2 and phase == "dc":
                phase = "newton"      # exchange is stuck; polish the reference
                stall = 0
            elif stall > 2:
                break
        else:
            stall = 0
        ref = new_ref
        if phase == "dc" and it < dc_iters:
            try:
                P, Q, E = dc_solve(m, k, ref, prob)
            except (ZeroDivisionError, ValueError) as ex:
                newton_note = f"dc singular: {ex!r}"
                break
        else:
            phase = "newton"
            if not newton_used:
                newton_used = True
            P2, Q2, E2, ok, note, res = newton_solve(
                m, k, ref, prob, P, Q, E)
            newton_ok = ok
            newton_note = note
            newton_res = res
            if ok:
                P, Q, E = P2, Q2, E2
            else:
                try:
                    P, Q, E = dc_solve(m, k, ref, prob)
                    phase = "dc"
                except (ZeroDivisionError, ValueError) as ex:
                    newton_note = f"dc singular after Newton: {ex!r}"
                    break
        errs, qvals, gm, g, signs_ok, ratio = consider(
            P, Q, E, ref, it, phase)

    chosen = cert if cert is not None else best
    if chosen is None:
        return {
            "m": m, "k": k, "stored": m + k + 1, "converged": False,
            "status": "no-admissible-fit", "certificate": False,
            "seconds": round(time.perf_counter() - started, 2),
            "newton_used": newton_used, "newton_ok": newton_ok,
            "newton_note": newton_note,
            "newton_residual": (mp.nstr(newton_res, 5)
                                if newton_res is not None else None),
            "history": [{key: mp.nstr(val, 8) if isinstance(val, mp.mpf)
                         else val for key, val in h.items()} for h in history],
        }
    entry, P, Q, E, ref, errs, qvals = chosen
    converged = cert is not None
    return {
        "m": m, "k": k, "stored": m + k + 1,
        "converged": converged,
        "status": ("converged" if converged else "best-admissible"),
        "certificate": converged,
        "grid_max": mp.nstr(entry["grid_max"], 17),
        "guard": mp.nstr(entry["guard"], 17),
        "ratio": mp.nstr(entry["ratio"], 17),
        "signs_ok": entry["signs_ok"],
        "E": mp.nstr(entry["E"], 17),
        "ref_min": mp.nstr(entry["ref_min"], 17),
        "ref": [mp.nstr(t, 25) for t in ref],
        "iterations": entry["it"],
        "seconds": round(time.perf_counter() - started, 2),
        "newton_used": newton_used, "newton_ok": newton_ok,
        "newton_note": newton_note,
        "newton_residual": (mp.nstr(newton_res, 5)
                            if newton_res is not None else None),
        "P": [mp.nstr(c, 70) for c in P],
        "Q": [mp.nstr(c, 70) for c in Q],
        "meets_raw": bool(entry["grid_max"] <= target
                          and entry["guard"] >= DEN_FLOOR),
        "history": [{key: mp.nstr(val, 8) if isinstance(val, mp.mpf)
                     else val for key, val in h.items()} for h in history],
    }


# ---------------------------------------------------------------------------
# 60-digit scoring, exactly the control's protocol


def score60(P, Q, a, W, target, n=2000, grid=None):
    if grid is None:
        grid = build_grid(a, W, n)
    ts, tsm, fs = grid
    with mp.workdps(EVAL_DPS):
        P60s = [mp.nstr(c, EVAL_DPS) for c in P]
        Q60s = [mp.nstr(c, EVAL_DPS) for c in Q]
        P60 = [mp.mpf(c) for c in P60s]
        Q60 = [mp.mpf(c) for c in Q60s]
        wr = mp.mpf(0)
        wa = mp.mpf(0)
        at = None
        dmin = mp.inf
        dmax = mp.mpf(0)
        for t, tmm, f in zip(ts, tsm, fs):
            fm = mp.mpf(f)
            num = poly_eval(P60, tmm)
            den = poly_eval(Q60, tmm)
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
    guard = dmin / dmax
    return {
        "e60": mp.nstr(wr, 17),
        "e60_abs": mp.nstr(wa, 17),
        "e60_at_t": repr(at),
        "guard60": mp.nstr(guard, 17),
        "meets60": bool(wr <= target and guard >= DEN_FLOOR),
        "P60": P60s,
        "Q60": Q60s,
    }


# ---------------------------------------------------------------------------
# Known-answer test: minimax polynomial for t^n on [-1,1]


def test_polynomial():
    print("known-answer test: best degree n-1 polynomial for t^n on [-1,1]")
    ok_all = True
    n = 2000
    ts = [2.0 * i / n - 1.0 for i in range(n + 1)]
    tsm = [mp.mpf(t) for t in ts]
    for n_deg in (4, 5, 6, 8):
        m = n_deg - 1
        prob = problem_power(n_deg)
        fs = [prob.f(t) for t in tsm]
        res = solve_pair(m, 0, tsm, fs, prob, mp.mpf("1e-300"))
        expected = mp.mpf(2) ** (1 - n_deg)
        got = mp.mpf(res["grid_max"])
        ok = (res["converged"]
              and abs(got - expected) <= mp.mpf("1e-12") * expected)
        coef_ok = None
        if n_deg == 4:
            P = [mp.mpf(c) for c in res["P"]]
            want = [mp.mpf("-0.125"), mp.mpf(0), mp.mpf(1), mp.mpf(0)]
            coef_ok = all(abs(P[j] - want[j]) <= mp.mpf("1e-40")
                          for j in range(4))
        ok_all = ok_all and ok and (coef_ok is not False)
        print(f"  t^{n_deg}: converged={res['converged']} "
              f"e={res['grid_max']} expected={mp.nstr(expected, 17)} "
              f"ratio={res['ratio']} coeff_ok={coef_ok} "
              f"-> {'OK' if ok else 'FAIL'}")
    print("test:", "PASS" if ok_all else "FAIL")
    return ok_all


# ---------------------------------------------------------------------------
# CLI


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--test", action="store_true")
    ap.add_argument("--a", default="62.0")
    ap.add_argument("--width", default=None)
    ap.add_argument("--m", type=int, default=None)
    ap.add_argument("--k", type=int, default=None)
    ap.add_argument("--out", default=None)
    ap.add_argument("--verbose", action="store_true")
    args = ap.parse_args(argv)
    if args.test:
        return 0 if test_polynomial() else 1
    if None in (args.width, args.m, args.k):
        raise SystemExit("need --width --m --k (or --test)")
    a = float(args.a)
    W = float(args.width)
    ts, tsm, fs = build_grid(a, W, 2000)
    prob = problem_f0(a, W)
    res = solve_pair(args.m, args.k, tsm, fs, prob, mp.mpf("1e-15"),
                     verbose=args.verbose)
    if "P" in res:
        P = [mp.mpf(c) for c in res["P"]]
        Q = [mp.mpf(c) for c in res["Q"]]
        res.update(score60(P, Q, a, W, 1e-15, grid=(ts, tsm, fs)))
        res.pop("P60", None)
        res.pop("Q60", None)
    out = json.dumps(res, indent=1)
    print(out)
    if args.out:
        Path(args.out).write_text(out, encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
