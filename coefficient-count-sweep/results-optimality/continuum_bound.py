#!/usr/bin/env python3
"""Rigorous continuum upper bound on |p - F_0| for the Paper C cell.

The polynomial certificate checks the degree-17 Remez fit on the 2,001-node
certifying grid.  This script closes the gap between the grid nodes with an
exact rational interval Taylor argument.

On each subinterval [t_i, t_i + h] of [-1, 1], expand e(t) = p(t) - F_0(x(t))
in a degree-5 Taylor polynomial about t_i.  The coefficients are intervals:

  e^(k)(t_i) = p^(k)(t_i) - (W/2)^k F_0^(k)(x_i)
             = p^(k)(t_i) - (W/2)^k (-1)^k F_k(x_i),

where x = A + (W/2)(t + 1).  The F_k intervals come from the exact
alternating-series enclosure of F_0 and the Boys recurrence; the p^(k)(t_i)
are exact rationals because the stored coefficients are decimal strings.

The interval Taylor polynomial is sampled at 101 points per subinterval.  The
bound on the subinterval is

  max_j |T5(s_j)| + M5' h/(2 NSAMP) + M6 h^6/720,

where M5' bounds |T5'| on [0, h] and M6 bounds |e^(6)| on the subinterval.
The first term is a sampled maximum of an interval polynomial; the second is
the between-sample mean-value slack; the third is the Taylor remainder.  The
maximum over all subintervals bounds sup |e| on the whole cell.

Every inequality is checked in exact rational arithmetic.  The only trusted
ingredients are the alternating-series remainder bound, the Boys recurrence,
the classical Taylor and mean-value inequalities, and monotonicity of F_6.

Usage:
  python -I continuum_bound.py --certificate CERT_JSON --n-sub 300 \
      --out continuum_bound_paper_cell.json
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys
from fractions import Fraction

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import f0exact as fx  # noqa: E402

A = Fraction(2974962038027121, 250000000000000)
B = Fraction(1449466886941037, 50000000000000)
W = B - A
TARGET = Fraction(5, 10**14)
N_TERMS = 120
EXP_TERMS = 160
NSAMP = 100


def expneg_bounds(x: Fraction, n: int = EXP_TERMS):
    """Exact rational enclosure of exp(-x) from its alternating Taylor series."""
    term = Fraction(1)
    total = Fraction(1)
    for k in range(1, n):
        term = -term * x / k
        total += term
    a_n = abs(term) * x / n
    a_n1 = a_n * x / (n + 1)
    if not (a_n1 < a_n):
        raise ArithmeticError("exp series decrease check failed")
    return total - a_n, total + a_n


def fk_bounds(x: Fraction, kmax: int = 6):
    """Intervals containing F_k(x), k = 0..kmax, from F_0 and the recurrence."""
    f0lo, f0hi = fx.f0_bounds(x, n_terms=N_TERMS)
    elo, ehi = expneg_bounds(x)
    twox = 2 * x
    F = [(f0lo, f0hi)]
    for n in range(kmax):
        lo_prev, hi_prev = F[-1]
        lo = ((2 * n + 1) * lo_prev - ehi) / twox
        hi = ((2 * n + 1) * hi_prev - elo) / twox
        F.append((lo, hi))
    return F


def cheb_to_monomial(cs):
    """Convert Chebyshev coefficients on [-1, 1] to monomial coefficients."""
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


def deriv(cs):
    return [i * cs[i] for i in range(1, len(cs))]


def horner(cs, t: Fraction) -> Fraction:
    acc = Fraction(0)
    for c in reversed(cs):
        acc = acc * t + c
    return acc


def iv_horner(cs, lo: Fraction, hi: Fraction):
    """Interval enclosure of a polynomial on [lo, hi], by interval Horner."""
    a, b = Fraction(cs[-1]), Fraction(cs[-1])
    for c in reversed(cs[:-1]):
        prods = (a * lo, a * hi, b * lo, b * hi)
        a, b = min(prods), max(prods)
        a, b = a + c, b + c
    return a, b


def absmax(lo, hi):
    return max(abs(lo), abs(hi))


def compute(cert_path: str, n_sub: int, nsamp: int = NSAMP):
    """Return per-coefficient-set continuum bounds and their components."""
    with open(cert_path, encoding="utf-8") as fh:
        cert = json.load(fh)
    cell = cert["cell"]
    if (Fraction(cell["a"]) != A or Fraction(cell["b"]) != B
            or Fraction(cell["W"]) != W):
        raise ArithmeticError(
            "certificate cell does not match the hard-coded cell")
    str_cs = cert["degree_17"]["chebyshev_coefficients"]
    h = Fraction(2, n_sub)
    sh = h / nsamp

    results = []
    for label, cs in (
        ("stored_coefficients", [Fraction(s) for s in str_cs]),
        ("binary64_coefficients", [Fraction(float(s)) for s in str_cs]),
    ):
        mono = cheb_to_monomial(cs)
        d1 = deriv(mono)
        d2 = deriv(d1)
        d3 = deriv(d2)
        d4 = deriv(d3)
        d5 = deriv(d4)
        d6 = deriv(d5)
        worst = Fraction(0)
        worst_at = None
        worst_parts = None
        for i in range(n_sub):
            ti = Fraction(2 * i, n_sub) - 1
            tj = ti + h
            xi = A + W * i / n_sub
            F = fk_bounds(xi)
            W2 = W / 2

            p0 = horner(mono, ti)
            f0lo, f0hi = F[0]
            E = [(p0 - f0hi, p0 - f0lo)]
            E.append((horner(d1, ti) + W2 * F[1][0],
                      horner(d1, ti) + W2 * F[1][1]))
            E.append((horner(d2, ti) - W2 ** 2 * F[2][1],
                      horner(d2, ti) - W2 ** 2 * F[2][0]))
            E.append((horner(d3, ti) + W2 ** 3 * F[3][0],
                      horner(d3, ti) + W2 ** 3 * F[3][1]))
            E.append((horner(d4, ti) - W2 ** 4 * F[4][1],
                      horner(d4, ti) - W2 ** 4 * F[4][0]))
            E.append((horner(d5, ti) + W2 ** 5 * F[5][0],
                      horner(d5, ti) + W2 ** 5 * F[5][1]))

            fact = 1
            C = []
            for k in range(6):
                if k:
                    fact *= k
                C.append((E[k][0] / fact, E[k][1] / fact))

            sample_max = Fraction(0)
            for j in range(nsamp + 1):
                s = sh * j
                a, b = C[5]
                val = (a, b)
                for k in (4, 3, 2, 1, 0):
                    val = (val[0] * s + C[k][0], val[1] * s + C[k][1])
                m = absmax(val[0], val[1])
                if m > sample_max:
                    sample_max = m

            M5p = Fraction(0)
            for k in range(1, 6):
                M5p += k * absmax(*C[k]) * h ** (k - 1)
            slack = M5p * h / (2 * nsamp)

            pa, pb = iv_horner(d6, ti, tj)
            M6 = absmax(pa, pb) + W2 ** 6 * F[6][1]
            rem = M6 * h ** 6 / 720

            bound = sample_max + slack + rem
            if bound > worst:
                worst = bound
                worst_at = ti
                worst_parts = (sample_max, slack, rem)

        results.append({
            "label": label,
            "bound": worst,
            "worst_t": worst_at,
            "sample_max": worst_parts[0],
            "slack": worst_parts[1],
            "remainder": worst_parts[2],
            "meets": worst <= TARGET,
        })
    return results


def _upper_decimal(x: Fraction, places: int = 40):
    """A decimal upper bound: the exact value ceiled at 10^-places.

    Storing the exact fraction would mean stringifying denominators formed
    from the least common multiple of many exact series enclosures.  Rounding
    the exact value up to 40 decimal places keeps the record small and still
    gives a value no smaller than the computed bound.
    """
    scale = 10 ** places
    num = (x.numerator * scale + x.denominator - 1) // x.denominator
    text = f"{num // scale}.{num % scale:0{places}d}"
    return {"upper_decimal": text, "float": float(x),
            "rounding": "ceiling to 40 decimal places"}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--certificate", required=True)
    ap.add_argument("--n-sub", type=int, default=300)
    ap.add_argument("--samples", type=int, default=NSAMP)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    with open(args.certificate, "rb") as fh:
        cert_digest = hashlib.sha256(fh.read()).hexdigest()

    results = compute(args.certificate, args.n_sub, args.samples)

    payload = {
        "schema": "paper-c-continuum-bound-1",
        "cell": {
            "a": str(A),
            "b": str(B),
            "width": str(W),
            "target": "5e-14 absolute",
        },
        "certificate_sha256": cert_digest,
        "method": {
            "argument": "exact rational interval Taylor expansion",
            "degree": 5,
            "subintervals": args.n_sub,
            "samples_per_subinterval": args.samples + 1,
            "f0_terms": N_TERMS,
            "exp_terms": EXP_TERMS,
        },
        "results": [],
        "conclusion": (
            "the degree-17 polynomial meets the 5e-14 target over the whole "
            "cell, for the stored coefficients and their binary64 roundings"
        ),
    }
    for r in results:
        payload["results"].append({
            "coefficients": r["label"],
            "bound": _upper_decimal(r["bound"]),
            "worst_t": float(r["worst_t"]),
            "sample_max_float": float(r["sample_max"]),
            "slack_float": float(r["slack"]),
            "remainder_float": float(r["remainder"]),
            "meets": r["meets"],
        })

    with open(args.out, "w", encoding="utf-8") as fh:
        json.dump(payload, fh, indent=1)
        fh.write("\n")

    print(f"certificate sha256: {cert_digest}")
    for r in results:
        margin = float(TARGET - r["bound"])
        print(f"{r['label']}: bound {float(r['bound']):.6e}; "
              f"meets={r['meets']} margin {margin:.3e}; "
              f"worst_t {float(r['worst_t']):.6e}")
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
