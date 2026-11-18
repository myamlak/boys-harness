#!/usr/bin/env python3
"""Secondary pilot: is the Chebyshev count 19 a true polynomial minimum?

At the paper cell (F_0 on [x0, x1], 5e-14 absolute) the paper's Chebyshev
pipeline stores 19 coefficients (degree 18).  The paper calls this a
solver-relative count.  This script asks the stronger polynomial-family
question:

  * upper bound: does some degree-17 polynomial meet 5e-14 on the 2,001-node
    grid (stored count 18)?
  * lower bound: can a de la Vallee Poussin certificate with a degree-16
    reference and 18 alternating extrema prove that no degree-16 polynomial
    meets?

A degree-d polynomial stores d+1 coefficients in either the monomial or the
Chebyshev basis, so the count question is degree-only.  The Remez exchange
below is a numerical construction; every reported inequality is rechecked in
exact rational arithmetic on polynomials with rational coefficients.

Usage: python -I certify_polynomial.py [--d17-only]
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time
from fractions import Fraction

import mpmath as mp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import f0exact as fx  # noqa: E402

A = Fraction("11.899848152108484")
B = Fraction("28.989337738820740")
W = B - A
TARGET = fx.TARGET
N_TERMS = 120
MP_DPS = 80


def f0_mp(x):
    return mp.sqrt(mp.pi) / (2 * mp.sqrt(x)) * mp.erf(mp.sqrt(x))


def f1_mp(x):
    return (f0_mp(x) - mp.e ** (-x)) / (2 * x)


def cheb_T(cs, t):
    """Evaluate a Chebyshev series by the three-term recurrence."""
    if not cs:
        return mp.mpf(0)
    b1 = mp.mpf(0)
    b2 = mp.mpf(0)
    for c in reversed(cs[1:]):
        b0 = 2 * t * b1 - b2 + c
        b2, b1 = b1, b0
    return t * b1 - b2 + cs[0]


def cheb_Tprime(cs, t):
    """Derivative of a Chebyshev series, via the differentiated recurrence."""
    d = len(cs) - 1
    if d <= 0:
        return mp.mpf(0)
    T = [mp.mpf(1), t]
    dT = [mp.mpf(0), mp.mpf(1)]
    for j in range(1, d):
        T.append(2 * t * T[j] - T[j - 1])
        dT.append(2 * T[j] + 2 * t * dT[j] - dT[j - 1])
    return sum(cs[j] * dT[j] for j in range(1, d + 1))


def cheb_to_monomial_exact(cs):
    """Exact conversion of Chebyshev coefficients to monomial coefficients."""
    d = len(cs) - 1
    T = {0: [Fraction(1)]}
    if d >= 1:
        T[1] = [Fraction(0), Fraction(1)]
    for j in range(1, d):
        prev, cur = T[j - 1], T[j]
        nxt = [Fraction(0)] * (j + 2)
        for i, c in enumerate(cur):
            nxt[i + 1] += 2 * c
        for i, c in enumerate(prev):
            nxt[i] -= c
        T[j + 1] = nxt
    out = [Fraction(0)] * (d + 1)
    for j, c in enumerate(cs):
        if c == 0:
            continue
        for i, tc in enumerate(T[j]):
            out[i] += c * tc
    return out


def cheb_eval_exact(cs, t):
    """Exact Chebyshev-series evaluation at a rational t."""
    if not cs:
        return Fraction(0)
    T0, T1 = Fraction(1), t
    total = cs[0]
    if len(cs) == 1:
        return total
    total += cs[1] * T1
    for j in range(1, len(cs) - 1):
        T2 = 2 * t * T1 - T0
        total += cs[j + 1] * T2
        T0, T1 = T1, T2
    return total


def remez_polynomial(d, n_iter=60, n_sample=8001):
    """Polynomial minimax on t in [-1,1] by the Chebyshev-basis Remez."""
    m = d + 2
    nodes = [-mp.cos(mp.pi * i / (d + 1)) for i in range(m)]

    def f(t):
        return f0_mp(mp.mpf(A.numerator) / A.denominator +
                     mp.mpf(W.numerator) / W.denominator * (t + 1) / 2)

    best = None
    for it in range(n_iter):
        A_mat = []
        rhs = []
        for i, t in enumerate(nodes):
            row = []
            T0, T1 = mp.mpf(1), t
            row.append(T0)
            if d >= 1:
                row.append(T1)
            for j in range(1, d):
                T2 = 2 * t * T1 - T0
                row.append(T2)
                T0, T1 = T1, T2
            row.append(mp.mpf((-1) ** i))
            A_mat.append(row)
            rhs.append(f(t))
        try:
            sol = mp.lu_solve(A_mat, rhs)
        except ZeroDivisionError:
            break
        cs = list(sol[:d + 1])
        err_nodes = abs(sol[d + 1])
        # sample the error and its derivative
        ts = [mp.mpf(2 * i) / (n_sample - 1) - 1 for i in range(n_sample)]
        vals = [cheb_T(cs, t) - f(t) for t in ts]
        dvals = []
        for t in ts:
            x = mp.mpf(A.numerator) / A.denominator + \
                mp.mpf(W.numerator) / W.denominator * (t + 1) / 2
            dvals.append(cheb_Tprime(cs, t) +
                         (mp.mpf(W.numerator) / W.denominator / 2) * f1_mp(x))
        pts = [(ts[0], vals[0]), (ts[-1], vals[-1])]
        for i in range(n_sample - 1):
            if dvals[i] == 0:
                pts.append((ts[i], vals[i]))
            elif (dvals[i] > 0) != (dvals[i + 1] > 0):
                lo, hi = ts[i], ts[i + 1]
                for _ in range(80):
                    mid = (lo + hi) / 2
                    xmid = mp.mpf(A.numerator) / A.denominator + \
                        mp.mpf(W.numerator) / W.denominator * (mid + 1) / 2
                    dmid = cheb_Tprime(cs, mid) + \
                        (mp.mpf(W.numerator) / W.denominator / 2) * f1_mp(xmid)
                    if (dmid > 0) == (dvals[i] > 0):
                        lo = mid
                    else:
                        hi = mid
                mid = (lo + hi) / 2
                pts.append((mid, cheb_T(cs, mid) - f(mid)))
        pts.sort(key=lambda z: z[0])
        # maximal alternating subsequence
        alt = []
        for t, v in pts:
            if not alt or (v > 0) != (alt[-1][1] > 0):
                alt.append((t, v))
            elif abs(v) > abs(alt[-1][1]):
                alt[-1] = (t, v)
        if len(alt) < m:
            break
        # contiguous window of length m with the largest minimum magnitude
        best_win, best_min = None, None
        for start in range(len(alt) - m + 1):
            win = alt[start:start + m]
            mn = min(abs(v) for _, v in win)
            if best_min is None or mn > best_min:
                best_win, best_min = win, mn
        new_nodes = [t for t, _ in best_win]
        shift = max(abs(new_nodes[i] - nodes[i]) for i in range(m))
        nodes = new_nodes
        max_err = max(abs(v) for _, v in pts)
        best = dict(d=d, cs=cs, E=err_nodes, max_err=max_err, it=it)
        if shift < mp.mpf("1e-13"):
            break
    return best


def exact_grid_poly(cs, n_nodes=2000, n_terms=N_TERMS, binary64=False):
    """Exact max |p(t)-F_0(x)| over the uniform t-grid."""
    if binary64:
        cs = [Fraction(float(c)) for c in cs]
    else:
        cs = [Fraction(mp.nstr(c, 60)) for c in cs]
    worst_hi, worst_lo = Fraction(0), Fraction(0)
    worst_t = None
    for i in range(n_nodes + 1):
        t = Fraction(2 * i, n_nodes) - 1
        x = A + W * (t + 1) / 2
        lo, hi = fx.f0_bounds(x, n_terms=n_terms)
        p = cheb_eval_exact(cs, t)
        e_hi = max(abs(lo - p), abs(hi - p))
        e_lo = Fraction(0) if lo <= p <= hi else min(abs(lo - p), abs(hi - p))
        if e_hi > worst_hi:
            worst_hi, worst_t = e_hi, t
        worst_lo = max(worst_lo, e_lo)
    return dict(max_error_upper=float(worst_hi),
                max_error_lower=float(worst_lo),
                worst_t=float(worst_t),
                meets=bool(worst_hi <= TARGET),
                target="5e-14")


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--d17-only", action="store_true")
    ap.add_argument("--out", default="certificate_polynomial_paper_cell.json")
    args = ap.parse_args(argv)
    started = time.time()
    mp.mp.dps = MP_DPS
    out = dict(
        schema="paper-c-polynomial-minimum-pilot-1",
        cell=dict(a=str(A), b=str(B), W=str(W), target="5e-14 absolute",
                  nodes=2001),
        note=("The paper's Chebyshev count is a Chebyshev-interpolation "
              "count. This probe answers the stronger question over all "
              "polynomials. A degree-17 polynomial stores 18 coefficients."),
    )
    ok17 = ok16 = None
    for d in ((17,) if args.d17_only else (17, 16)):
        t0 = time.time()
        rem = remez_polynomial(d)
        if rem is None:
            print(f"degree {d}: Remez failed")
            out[f"degree_{d}"] = dict(remez=False)
            continue
        print(f"degree {d}: Remez iter={rem['it']} E={mp.nstr(rem['E'], 6)} "
              f"max={mp.nstr(rem['max_err'], 6)} "
              f"({time.time()-t0:.1f}s)")
        if d == 17:
            hp = exact_grid_poly(rem["cs"])
            b64 = exact_grid_poly(rem["cs"], binary64=True)
            print(f"  exact grid: hp={hp['max_error_upper']:.6e} "
                  f"meets={hp['meets']}; binary64 coeffs="
                  f"{b64['max_error_upper']:.6e} meets={b64['meets']}")
            out["degree_17"] = dict(
                remez_E=float(rem["E"]), remez_max=float(rem["max_err"]),
                chebyshev_coefficients=[mp.nstr(c, 40) for c in rem["cs"]],
                exact_grid_high_precision=hp,
                exact_grid_binary64_coefficients=b64,
                stored_count=18,
                conclusion=("a degree-17 polynomial meets 5e-14"
                            if hp["meets"] else
                            "the degree-17 Remez fit misses 5e-14"),
            )
            ok17 = hp["meets"]
        else:
            # Exact de la Vallee Poussin certificate for degree 16:
            # 18 alternating extrema, checked with exact rational polynomials.
            cs = [Fraction(mp.nstr(c, 50)) for c in rem["cs"]]
            # Use the Remez extrema as nodes: reconstruct them by sampling.
            nodes = []
            m = d + 2
            cs_mp = [mp.mpf(mp.nstr(c, 50)) for c in rem["cs"]]
            # reuse the alternating extrema found in the last Remez iterate
            # by re-running the sampler once
            ts = [mp.mpf(2 * i) / 8000 - 1 for i in range(8001)]
            def f(t):
                return f0_mp(mp.mpf(A.numerator) / A.denominator +
                             mp.mpf(W.numerator) / W.denominator *
                             (t + 1) / 2)
            vals = [cheb_T(cs_mp, t) - f(t) for t in ts]
            dvals = []
            for t in ts:
                x = mp.mpf(A.numerator) / A.denominator + \
                    mp.mpf(W.numerator) / W.denominator * (t + 1) / 2
                dvals.append(cheb_Tprime(cs_mp, t) +
                             (mp.mpf(W.numerator) / W.denominator / 2) *
                             f1_mp(x))
            pts = [(ts[0], vals[0]), (ts[-1], vals[-1])]
            for i in range(len(ts) - 1):
                if dvals[i] == 0:
                    pts.append((ts[i], vals[i]))
                elif (dvals[i] > 0) != (dvals[i + 1] > 0):
                    lo, hi = ts[i], ts[i + 1]
                    for _ in range(90):
                        mid = (lo + hi) / 2
                        xmid = mp.mpf(A.numerator) / A.denominator + \
                            mp.mpf(W.numerator) / W.denominator * \
                            (mid + 1) / 2
                        dmid = cheb_Tprime(cs_mp, mid) + \
                            (mp.mpf(W.numerator) / W.denominator / 2) * \
                            f1_mp(xmid)
                        if (dmid > 0) == (dvals[i] > 0):
                            lo = mid
                        else:
                            hi = mid
                    mid = (lo + hi) / 2
                    pts.append((mid, cheb_T(cs_mp, mid) - f(mid)))
            pts.sort(key=lambda z: z[0])
            alt = []
            for t, v in pts:
                if not alt or (v > 0) != (alt[-1][1] > 0):
                    alt.append((t, v))
                elif abs(v) > abs(alt[-1][1]):
                    alt[-1] = (t, v)
            if len(alt) < m:
                print(f"  dVP degree 16: only {len(alt)} alternating "
                      f"extrema, need {m}")
                out["degree_16"] = dict(remez_E=float(rem["E"]),
                                        remez_max=float(rem["max_err"]),
                                        certified=False,
                                        error="too few alternations")
                ok16 = False
                continue
            # choose a contiguous window of length m
            best_win, best_min = None, None
            for start in range(len(alt) - m + 1):
                win = alt[start:start + m]
                mn = min(abs(v) for _, v in win)
                if best_min is None or mn > best_min:
                    best_win, best_min = win, mn
            for t, _v in best_win:
                x = mp.mpf(A.numerator) / A.denominator + \
                    mp.mpf(W.numerator) / W.denominator * (t + 1) / 2
                xf = Fraction(mp.nstr(x, 20))
                tf = 2 * (xf - A) / W - 1
                nodes.append((xf, tf))
            P = cheb_to_monomial_exact(cs)
            Q = [Fraction(1)]
            okd, delta, rows = fx.dvp_check_mapped(P, Q, nodes,
                                                   target=TARGET,
                                                   n_terms=N_TERMS)
            print(f"  dVP degree 16: ok={okd} delta={float(delta):.6e} "
                  f"({float(delta)/float(TARGET):.2f} x target) "
                  f"nodes={len(nodes)}")
            out["degree_16"] = dict(
                remez_E=float(rem["E"]), remez_max=float(rem["max_err"]),
                needed_alternations=m, nodes=rows,
                certified=bool(okd),
                delta="%.6e" % float(delta),
                delta_over_target="%.6f" % (float(delta)/float(TARGET)),
                P_monomial=["%d/%d" % (c.numerator, c.denominator)
                            for c in P],
                Q=[str(c) for c in Q],
            )
            ok16 = okd
    out["conclusion"] = dict(
        polynomial_minimum_stored_count=(
            18 if (ok17 and ok16) else None),
        polynomial_minimum_upper_bound=18 if ok17 else None,
        degree_17_meets=bool(ok17),
        degree_16_certified_to_fail=bool(ok16),
        paper_chebyshev_count=19,
        note=("The paper's 19 is a Chebyshev-interpolation count, not a "
              "polynomial-family minimum; the paper labels it solver-relative."
              if (ok17 and ok16) else
              "Incomplete: see the individual fields."),
    )
    with open(args.out, "w", encoding="utf-8") as fh:
        json.dump(out, fh, indent=1)
    print(f"wrote {args.out}")
    print(json.dumps(out["conclusion"], indent=1))
    return 0 if (ok17 and (ok16 or args.d17_only)) else 1


if __name__ == "__main__":
    raise SystemExit(main())
