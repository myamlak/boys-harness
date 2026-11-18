#!/usr/bin/env python3
"""Synthetic timing probe for a wide-cell dVP certificate (NOT a certificate).

The paper cell's count-11 types have total degree 10, so the pilot's
per-type cost is small.  This probe measures the same machinery on a type
with total degree 27 (stored count 28), the regime of the widest cells in the
10^-14 arm, where the count-below types are [27/1] ... [15/12].

It reports the BRASIL solve time, the monomial conversion, the numerical
extrema search, and the exact verification time, or the failure mode.  No
number here is a proof; the output is a cost estimate.

Usage: python -I bench_high_degree.py
"""

from __future__ import annotations

import os
import sys
import time
from fractions import Fraction

import mpmath as mp
import numpy as np
from baryrat import brasil
from scipy.special import erf

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import f0exact as fx  # noqa: E402
import make_certificate as mc  # noqa: E402

A = Fraction("11.899848152108484")
W = Fraction("19.100151847891516")  # the ship width, a wide cell
A0 = A
mc.A = A
mc.B = A + W
mc.W_RECORD = W
mc.W_PAPER = W


def f0_np(x):
    x = np.asarray(x, dtype=float)
    return np.sqrt(np.pi) / (2.0 * np.sqrt(x)) * erf(np.sqrt(x))


def main(m=15, k=12):
    mp.mp.dps = 60
    print(f"synthetic wide-cell probe: [{m}/{k}] stored {m+1+k}, "
          f"W={W} (ship width)")
    f = lambda t: f0_np(float(A) + float(W) *
                       (np.asarray(t) + 1.0) / 2.0)
    t0 = time.perf_counter()
    try:
        r, info = brasil(f, (-1.0, 1.0), (m, k), info=True,
                         tol=1e-4, maxiter=1000, npi=-30,
                         init_steps=100, max_step_size=0.1)
    except Exception as exc:  # noqa: BLE001
        print(f"BRASIL raised {exc!r}")
        return 1
    t_solve = time.perf_counter() - t0
    print(f"BRASIL: converged={info.converged} dev={info.deviation:.3e} "
          f"{t_solve:.1f}s")
    conv = mc.bary_to_monomial(r)
    if conv is None:
        print("monomial conversion failed")
        return 1
    P, Q = conv
    P = [mp.nstr(c, 60) for c in P]
    Q = [mp.nstr(c, 60) for c in Q]
    t0 = time.perf_counter()
    cert = mc.certify_type(P, Q, m, k, f"synthetic [{m}/{k}]", "benchmark")
    t_cert = time.perf_counter() - t0
    print(f"certify_type: certified={cert['certified']} "
          f"pole_free={cert['pole_free']} "
          f"nodes={len(cert.get('nodes', []))} "
          f"delta={cert.get('delta')} {t_cert:.1f}s")
    return 0 if cert["certified"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
