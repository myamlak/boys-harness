#!/usr/bin/env python3
"""Independent DCT check of the Paper C Chebyshev-interpolation counts.

The paper reports a Chebyshev-interpolation count of 19 stored coefficients
(degree 18) for the paper cell.  This script records the independent DCT
observation behind the paper's solver-relative label:

  * the degree-17 Chebyshev-Gauss (DCT-II) interpolant of F_0 on the paper
    cell misses the 5e-14 absolute target;
  * the degree-18 interpolant meets it.

The interpolant of degree d is the unique polynomial through the d+1
Chebyshev-Gauss nodes t_k = cos(pi (2k+1) / (2(d+1))), k = 0..d, of the
mapped variable t = 2 (x - x0) / W - 1.  Its Chebyshev coefficients are the
DCT-II sums

    c_j = (2/n) sum_k F_0(x(t_k)) cos(pi j (2k+1) / (2n)),  n = d + 1,

with c_0 halved, and p(t) = sum_j c_j T_j(t).

Everything here is computed from F_0 through its closed form in mpmath at
60 decimal digits; no coefficient is copied from the paper or from the
certificate deposit.  The error is measured on the same 2,001-node uniform
x-grid as the certificate, using the same cell endpoints.  A second pass at
80 digits confirms that the reported digits are converged.

The reproduction targets come from the independent review: degree 17 =
5.5998e-14, degree 18 = 9.88e-15.  The script does not assume those values;
it computes each maximum, then checks agreement to the reviewer's stated
precision.

Usage, from the results-optimality directory:
    python -I dct/dct_interpolation_check.py
Exit 0 iff the degree-17 miss, the degree-18 meeting, the agreement with
the review values, and the precision-stability check all hold.
"""

from __future__ import annotations

import hashlib
import json
import os
import sys
from fractions import Fraction

import mpmath as mp

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
RESULTS_DIR = os.path.dirname(SCRIPT_DIR)
CERT_PATH = os.path.join(RESULTS_DIR, "certificate_rational_paper_cell.json")
OUT_PATH = os.path.join(SCRIPT_DIR, "dct_interpolation_paper_cell.json")

DPS = 60
CHECK_DPS = 80
N_GRID = 2001
COMMAND = "python -I dct/dct_interpolation_check.py"

# Independent review values, quoted at the precision the reviewer stated.
REVIEW = {
    17: (mp.mpf("5.5998e-14"), mp.mpf("1e-3")),
    18: (mp.mpf("9.88e-15"), mp.mpf("1e-3")),
}
REVIEW_DEGREE = {17: False, 18: True}


def load_cell():
    """Read the paper cell from the certificate deposit."""
    with open(CERT_PATH, encoding="utf-8") as fh:
        cert = json.load(fh)
    cell = cert["cell"]
    x0 = Fraction(cell["x0"])
    x1 = Fraction(cell["x1"])
    W = Fraction(cell["W_paper"])
    if x1 - x0 != W:
        raise SystemExit("certificate W_paper is not x1 - x0")
    target = Fraction(5, 10**14)
    if cell.get("certifying_nodes") != N_GRID:
        raise SystemExit("unexpected certifying node count")
    with open(CERT_PATH, "rb") as fh:
        digest = hashlib.sha256(fh.read()).hexdigest()
    return x0, x1, W, target, digest


def f0_mp(x):
    return mp.sqrt(mp.pi) / (2 * mp.sqrt(x)) * mp.erf(mp.sqrt(x))


def dct_coefficients(degree, x0, W):
    """Chebyshev-Gauss DCT-II coefficients of the degree-d interpolant."""
    n = degree + 1
    x0m = mp.mpf(x0.numerator) / x0.denominator
    Wm = mp.mpf(W.numerator) / W.denominator
    thetas = [mp.pi * (2 * k + 1) / (2 * n) for k in range(n)]
    fvals = [f0_mp(x0m + Wm * (mp.cos(th) + 1) / 2) for th in thetas]
    cs = []
    for j in range(n):
        s = mp.fsum(fvals[k] * mp.cos(j * thetas[k]) for k in range(n))
        c = 2 * s / n
        if j == 0:
            c /= 2
        cs.append(c)
    return cs


def cheb_eval(cs, t):
    """Clenshaw evaluation of a Chebyshev series."""
    b1 = mp.mpf(0)
    b2 = mp.mpf(0)
    for c in reversed(cs[1:]):
        b0 = 2 * t * b1 - b2 + c
        b2, b1 = b1, b0
    return t * b1 - b2 + cs[0]


def grid_error(cs, x0m, Wm):
    """Max |p(t) - F_0(x)| on the uniform 2,001-node paper grid."""
    worst = mp.mpf(0)
    worst_i = 0
    for i in range(N_GRID):
        t = mp.mpf(2 * i) / (N_GRID - 1) - 1
        x = x0m + Wm * (t + 1) / 2
        err = abs(cheb_eval(cs, t) - f0_mp(x))
        if err > worst:
            worst, worst_i = err, i
    t = mp.mpf(2 * worst_i) / (N_GRID - 1) - 1
    return worst, worst_i, t


def compute(degree, x0, W):
    cs = dct_coefficients(degree, x0, W)
    x0m = mp.mpf(x0.numerator) / x0.denominator
    Wm = mp.mpf(W.numerator) / W.denominator
    worst, idx, t = grid_error(cs, x0m, Wm)
    return cs, worst, idx, t


def main():
    x0, x1, W, target, cert_digest = load_cell()
    target_mp = mp.mpf(target.numerator) / target.denominator
    print(f"command: {COMMAND}")
    print(f"working precision: {DPS} decimal digits "
          f"(stability check at {CHECK_DPS})")
    print(f"python: {sys.version.split()[0]}, mpmath: {mp.__version__}")
    print(f"cell: F_0(x) = sqrt(pi)/(2 sqrt(x)) erf(sqrt(x)) on [x0, x1]")
    print(f"  x0 = {x0} = {float(x0)!r}")
    print(f"  x1 = {x1} = {float(x1)!r}")
    print(f"  W  = {W} = {float(W)!r}")
    print(f"  target = {target} = {float(target):.6e} absolute")
    print(f"  grid = {N_GRID} uniform x-nodes (same as the certificate)")
    print(f"cell source: certificate_rational_paper_cell.json, "
          f"sha256 {cert_digest}")
    print("coefficient source: DCT-II of F_0 at the Chebyshev-Gauss nodes, "
          "computed here; none copied")

    out = {
        "schema": "paper-c-dct-interpolation-check-1",
        "command": COMMAND,
        "versions": {
            "python": sys.version.split()[0],
            "mpmath": mp.__version__,
        },
        "precision": {
            "working_dps": DPS,
            "stability_check_dps": CHECK_DPS,
        },
        "cell_source": {
            "path": "../certificate_rational_paper_cell.json",
            "sha256": cert_digest,
        },
        "cell": {
            "function": "F_0(x) = sqrt(pi)/(2 sqrt(x)) erf(sqrt(x))",
            "x0": f"{x0.numerator}/{x0.denominator}",
            "x1": f"{x1.numerator}/{x1.denominator}",
            "W": f"{W.numerator}/{W.denominator}",
            "target": "5e-14 absolute",
            "grid_nodes": N_GRID,
            "grid": "uniform in x on [x0, x1], including both endpoints",
        },
        "method": {
            "nodes": "Chebyshev-Gauss (first kind): "
                     "t_k = cos(pi (2k+1) / (2(d+1))), k = 0..d",
            "coefficients": "DCT-II: c_j = (2/n) sum_k f_k cos(pi j "
                            "(2k+1)/(2n)), n = d+1, c_0 halved",
            "evaluation": "Clenshaw recurrence on the 2,001-node grid",
            "coefficient_source": "computed in this script from the F_0 "
                                  "closed form at the stated precision; not "
                                  "copied from the paper or a certificate",
        },
        "independent_review": {
            "degree_17": "5.5998e-14",
            "degree_18": "9.88e-15",
            "agreement_tolerance_relative": 1e-3,
        },
        "results": {},
        "conclusion": {},
    }
    ok = True
    mp.mp.dps = DPS
    for degree in (17, 18):
        cs, worst, idx, t = compute(degree, x0, W)
        meets = bool(worst <= target_mp)
        value = mp.nstr(worst, 25)
        # Stability across working precision.
        mp.mp.dps = CHECK_DPS
        cs_hi, worst_hi, idx_hi, t_hi = compute(degree, x0, W)
        delta = abs(worst - worst_hi)
        stable = bool(delta < mp.mpf("1e-45"))
        # Coefficient drift between the two precisions.
        coeff_drift = max(abs(a - b) for a, b in zip(cs, cs_hi))
        expected, rel_tol = REVIEW[degree]
        rel = abs(worst - expected) / expected
        review_ok = bool(rel <= rel_tol)
        want_meets = REVIEW_DEGREE[degree]
        statement_ok = bool(meets == want_meets)
        ok &= statement_ok and review_ok and stable
        x_worst = mp.mpf(x0.numerator) / x0.denominator + \
            mp.mpf(W.numerator) / W.denominator * (t + 1) / 2
        print(f"\ndegree {degree}: n_coeffs={degree + 1}, stored_count="
              f"{degree + 1}")
        print(f"  coefficients (60-digit): "
              + ", ".join(mp.nstr(c, 50) for c in cs))
        print(f"  max |err| = {value} at t={mp.nstr(t, 8)}, "
              f"x={mp.nstr(x_worst, 20)} (grid node {idx})")
        print(f"  meets target: {meets} (expected {want_meets})")
        print(f"  independent review {mp.nstr(expected, 8)}: relative "
              f"deviation {mp.nstr(rel, 6)} (tolerance {mp.nstr(rel_tol, 2)})"
              f" -> {review_ok}")
        print(f"  stability: |dps{DPS} - dps{CHECK_DPS}| = "
              f"{mp.nstr(delta, 6)}, coefficient drift "
              f"{mp.nstr(coeff_drift, 6)} -> {stable}")
        out["results"][f"degree_{degree}"] = {
            "degree": degree,
            "n_coefficients": degree + 1,
            "stored_count": degree + 1,
            "chebyshev_coefficients_dps60": [mp.nstr(c, 50) for c in cs],
            "max_abs_error": value,
            "max_abs_error_float": float(worst),
            "argmax_t": mp.nstr(t, 12),
            "argmax_x": mp.nstr(x_worst, 30),
            "argmax_grid_index": idx,
            "target": "5e-14 absolute",
            "meets_target": meets,
            "review_value": mp.nstr(expected, 8),
            "review_relative_deviation": mp.nstr(rel, 6),
            "review_matches": review_ok,
            "precision_stability": {
                "check_dps": CHECK_DPS,
                "abs_difference": mp.nstr(delta, 6),
                "coefficient_drift": mp.nstr(coeff_drift, 6),
                "stable": stable,
            },
        }
        mp.mp.dps = DPS
    out["conclusion"] = {
        "degree_17_misses_target": not out["results"]["degree_17"]
        ["meets_target"],
        "degree_18_meets_target": out["results"]["degree_18"]
        ["meets_target"],
        "reproduces_independent_review": all(
            r["review_matches"] for r in out["results"].values()),
        "note": "The degree-17 DCT interpolant misses the target and the "
                "degree-18 interpolant meets it, reproducing the independent "
                "review; the paper's 19 is therefore a Chebyshev-"
                "interpolation (solver-relative) count, not a "
                "polynomial-family minimum.",
        "pass": bool(ok),
    }
    with open(OUT_PATH, "w", encoding="utf-8", newline="\n") as fh:
        json.dump(out, fh, indent=1)
        fh.write("\n")
    rel = os.path.relpath(OUT_PATH, RESULTS_DIR).replace(os.sep, "/")
    print(f"\nwrote {rel}")
    print(f"{'PASS' if ok else 'FAIL'}")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
