#!/usr/bin/env python3
"""Independent re-check of the Paper C rational optimality certificate.

Reads certificate_rational_paper_cell.json, rebuilds every exact rational
from the stored decimal strings, and re-verifies:
  * the stored 12-coefficient nominal truncation meets 5e-14 on the
    2,001-node paper grid, for both the 70-digit coefficients and the
    binary64-stored coefficients;
  * every recorded stored-11 type has a pole-free reference and m+k+2
    exactly alternating nodes whose exact lower bounds exceed 5e-14;
  * the search-excluded endpoint types, if present, do the same;
  * the conclusion fields agree with those checks.

A numerical mpmath cross-check of the exact F_0 interval is included at every
node (mpmath is not trusted for a certified number; it is a bug check on the
exact evaluator).

Usage: python -I verify_certificate.py [certificate.json]
Exit 0 iff every check passes.
"""

from __future__ import annotations

import json
import os
import sys
import time
from fractions import Fraction

import mpmath as mp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import f0exact as fx  # noqa: E402

N_TERMS = 200  # tighter than the generator's 120
TARGET = fx.TARGET


def fail(msg):
    print(f"FAIL: {msg}")
    return False


def frac_list(strs):
    return [Fraction(s) for s in strs]


def exact_grid(P, Q, W, target, n_terms=N_TERMS):
    worst_hi, worst_lo = Fraction(0), Fraction(0)
    for i in range(2001):
        t = Fraction(2 * i, 2000) - 1
        x = Fraction("11.899848152108484") + W * (t + 1) / 2
        lo, hi = fx.f0_bounds(x, n_terms=n_terms)
        r, den = fx.rat_eval(P, Q, t)
        dlo, dhi = lo - r, hi - r
        e_hi = max(abs(dlo), abs(dhi))
        e_lo = Fraction(0) if dlo <= 0 <= dhi else min(abs(dlo), abs(dhi))
        worst_hi = max(worst_hi, e_hi)
        worst_lo = max(worst_lo, e_lo)
    return worst_lo, worst_hi, worst_hi <= target


def main(argv):
    path = argv[1] if len(argv) > 1 else "certificate_rational_paper_cell.json"
    print(f"verifying {path}")
    with open(path, encoding="utf-8") as fh:
        cert = json.load(fh)
    if cert.get("schema") != "paper-c-minimax-certificate-pilot-1":
        return fail("unexpected schema") or 1
    target = TARGET
    ok = True
    t_start = time.time()

    # --- stored candidate upper-bound witnesses ---------------------------
    nominal = cert["stored_candidate"]
    P_hp = frac_list(nominal["P"])
    Q_hp = frac_list(nominal["Q"])
    W_paper = Fraction(cert["cell"]["W_paper"])
    lo, hi, meets = exact_grid(P_hp, Q_hp, W_paper, target)
    print(f"stored candidate (70-digit coeffs): max|err| in "
          f"[{float(lo):.6e}, {float(hi):.6e}], meets={meets}")
    ok &= meets

    b64 = cert.get("stored_candidate_binary64")
    if b64:
        P64 = frac_list(b64["P"])
        Q64 = frac_list(b64["Q"])
        lo64, hi64, meets64 = exact_grid(P64, Q64, W_paper, target)
        print(f"stored candidate (binary64 coeffs): max|err| in "
              f"[{float(lo64):.6e}, {float(hi64):.6e}], meets={meets64}")
        ok &= meets64

    # --- de la Vallee Poussin certificates --------------------------------
    recorded = cert["recorded_below_types"]
    expected = {(10 - k, k) for k in range(1, 10)}
    got = {(t["m"], t["k"]) for t in recorded}
    if got != expected:
        ok = fail(f"recorded type set {sorted(got)} != {sorted(expected)}")
    mp.mp.dps = 80
    for group, label in ((recorded, "recorded"), (cert.get("extra_types", []),
                                                  "extra")):
        for t in group:
            m, k = t["m"], t["k"]
            P = frac_list(t["P"])
            Q = frac_list(t["Q"])
            needed = m + k + 2
            if t["stored"] != m + 1 + k:
                ok = fail(f"{t['label']}: stored count mismatch")
            if len(t["nodes"]) < needed:
                ok = fail(f"{t['label']}: {len(t['nodes'])} nodes < {needed}")
                continue
            # pole-free check
            lo_t = Fraction(-1)
            hi_t = 2 * (Fraction(cert["cell"]["x1"]) -
                        Fraction(cert["cell"]["x0"])) / \
                Fraction(cert["cell"]["W_record"]) - 1
            if fx.poly_eval(Q, lo_t) == 0 or fx.poly_eval(Q, hi_t) == 0:
                ok = fail(f"{t['label']}: reference pole at an endpoint")
                continue
            roots = fx.sturm_root_count(Q, lo_t, hi_t)
            if roots != 0:
                ok = fail(f"{t['label']}: {roots} reference poles in interval")
                continue
            # exact re-verification at every node
            signs = []
            delta = None
            for node in t["nodes"]:
                x = Fraction(int(node["x_num"]), int(node["x_den"]))
                tt = Fraction(int(node["t_num"]), int(node["t_den"]))
                flo, fhi = fx.f0_bounds(x, n_terms=N_TERMS)
                r, _den = fx.rat_eval(P, Q, tt)
                dlo, dhi = flo - r, fhi - r
                if dlo > 0 and dhi > 0:
                    sign, mag = 1, dlo
                elif dlo < 0 and dhi < 0:
                    sign, mag = -1, -dhi
                else:
                    ok = fail(f"{t['label']}: indeterminate sign at x={x}")
                    break
                if sign != node["sign"]:
                    ok = fail(f"{t['label']}: stored sign mismatch at x={x}")
                    break
                stored = Fraction(int(node["lower_num"]),
                                  int(node["lower_den"]))
                if stored > mag:
                    ok = fail(f"{t['label']}: stored lower bound exceeds "
                              f"recomputed exact bound at x={x}")
                    break
                if mag <= target:
                    ok = fail(f"{t['label']}: magnitude {float(mag):.3e} "
                              f"<= target at x={x}")
                    break
                signs.append(sign)
                delta = mag if delta is None else min(delta, mag)
                # mpmath cross-check of the exact interval (bug check only;
                # the exact interval can be far narrower than mpmath's
                # working-precision error, so allow 1e-60 slack).
                xm = mp.mpf(x.numerator) / x.denominator
                vm = mp.sqrt(mp.pi) / (2 * mp.sqrt(xm)) * mp.erf(mp.sqrt(xm))
                lom = mp.mpf(flo.numerator) / flo.denominator
                him = mp.mpf(fhi.numerator) / fhi.denominator
                slack = mp.mpf("1e-60")
                if not (lom - slack <= vm <= him + slack):
                    ok = fail(f"{t['label']}: mpmath cross-check at x={x}")
                    break
            alt = all(signs[i] != signs[i + 1] for i in range(len(signs) - 1))
            if not alt:
                ok = fail(f"{t['label']}: signs do not alternate")
            if delta is None:
                ok = fail(f"{t['label']}: no verified node")
                continue
            print(f"{label:8s} [{m}/{k}] nodes={len(t['nodes'])} "
                  f"needed={needed} alternating={alt} "
                  f"delta={float(delta):.6e} "
                  f"({float(delta)/float(target):.2f} x target)")
    # --- extras coverage ---------------------------------------------------
    if cert["conclusion"].get("beyond_search_certified"):
        extra_set = {(t["m"], t["k"]) for t in cert.get("extra_types", [])}
        if not {(10, 0), (0, 10)} <= extra_set:
            ok = fail("beyond_search_certified but endpoint extras missing")

    # --- conclusion consistency -------------------------------------------
    conc = cert["conclusion"]
    if conc.get("recorded_family_certified") and not all(
            t["certified"] for t in recorded):
        ok = fail("conclusion claims recorded family but a type is uncertified")
    if conc.get("rational_minimum_stored_count") == 12 and not (
            all(t["certified"] for t in recorded)):
        ok = fail("conclusion claims minimum 12 without all types certified")
    if cert.get("defect_list"):
        ok = fail(f"defect list non-empty: {cert['defect_list']}")

    elapsed = time.time() - t_start
    print(f"{'PASS' if ok else 'FAIL'} in {elapsed:.1f}s")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
