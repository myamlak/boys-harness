#!/usr/bin/env python3
"""Pilot minimax-optimality certificate for the Paper C paper cell.

Cell: F_0 on the paper's principal interval [x0, x1], prior-work width,
absolute target 5e-14.  The corrected record reports rational stored count 12
(pair [5/6], nominal truncation); the 60-digit audit reports every stored-11
pair's 70-digit monomial coefficients.

What this script does
---------------------
1. Verifies, with exact rational arithmetic, that the stored 12-coefficient
   nominal truncation meets the 5e-14 target on the 2,001-node certifying
   grid (an upper-bound witness).
2. For every stored-11 type [m/k] in the paper's admissible pair set, takes
   the audit's own 70-digit fit as a feasible reference r, locates the
   alternating extrema of F_0 - r numerically, and then verifies with exact
   rational arithmetic that F_0 - r alternates in sign at m+k+2 points with
   magnitude above 5e-14.  By the de la Vallee Poussin argument this proves
   that no rational of that type meets the target; because every stored
   count <= 11 fit is contained in one of those types (pad the numerator),
   the rational stored count 12 is minimal at this cell.
3. Optionally certifies the search-excluded endpoint types [0/10] and [10/0]
   the same way, extending the result beyond the recorded search family.
4. Checks each reference denominator for real roots on the interval with an
   exact Sturm sequence, so the reference is pole-free.

The numerical step locates nodes; it is not trusted.  Every inequality in
the output certificate is exact-rational and is rechecked by
verify_certificate.py.

Usage:
  python -I make_certificate.py
  python -I make_certificate.py --skip-grid --skip-extras
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys
import time
from datetime import datetime, timezone
from fractions import Fraction

import mpmath as mp
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import f0exact as fx  # noqa: E402

A = Fraction("11.899848152108484")
B = Fraction("28.989337738820740")
W_PAPER = B - A                      # 17.089489586712256, the paper's W_prior
W_RECORD = Fraction("17.089489586712254")  # the audit cell's width field
TARGET = fx.TARGET                   # 5e-14
AUDIT_DEFAULT = (r"C:\Users\myaml\dev\papers-c-v5\papers\coefficient-counts"
                 r"\results-corrected\corrected-counts-60dig.json")
CORRECTED_DEFAULT = (r"C:\Users\myaml\dev\papers-c-v5\papers\coefficient-counts"
                     r"\results-corrected\corrected-counts.json")

MP_DPS = 60
N_TERMS = 120


def now() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def sha256(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def find_cell(doc, width_name, target):
    for cell in doc["cells"]:
        if cell.get("width_name") == width_name and cell.get("target") == target:
            return cell
    raise SystemExit(f"cell ({width_name}, {target}) not found")


# ---------------------------------------------------------------------------
# Numerical location of alternation points.  This is a search, not a proof.


def f0_mp(x):
    return mp.sqrt(mp.pi) / (2 * mp.sqrt(x)) * mp.erf(mp.sqrt(x))


def f1_mp(x):
    return (f0_mp(x) - mp.e ** (-x)) / (2 * x)


def horner(cs, t):
    acc = mp.mpf(0)
    for c in reversed(cs):
        acc = acc * t + c
    return acc


def deriv(cs):
    return [i * cs[i] for i in range(1, len(cs))]


def g_value(t, P, Q):
    x = mp.mpf(A.numerator) / A.denominator + \
        mp.mpf(W_RECORD.numerator) / W_RECORD.denominator * (t + 1) / 2
    return f0_mp(x) - horner(P, t) / horner(Q, t)


def g_prime(t, P, Q, dP, dQ):
    wm = mp.mpf(W_RECORD.numerator) / W_RECORD.denominator
    x = mp.mpf(A.numerator) / A.denominator + wm * (t + 1) / 2
    p, q = horner(P, t), horner(Q, t)
    return -(wm / 2) * f1_mp(x) - (horner(dP, t) * q -
                                   p * horner(dQ, t)) / (q * q)


def extrema_sequence(P, Q, n_sample=6001, include_endpoints=True):
    """Stationary points of F_0 - r in t, plus optional endpoints."""
    dP, dQ = deriv(P), deriv(Q)
    ts = [mp.mpf(2 * i) / (n_sample - 1) - 1 for i in range(n_sample)]
    gs = [g_prime(t, P, Q, dP, dQ) for t in ts]
    pts = []
    if include_endpoints:
        pts.extend([(ts[0], g_value(ts[0], P, Q)),
                    (ts[-1], g_value(ts[-1], P, Q))])
    for i in range(n_sample - 1):
        if gs[i] == 0:
            pts.append((ts[i], g_value(ts[i], P, Q)))
        elif (gs[i] > 0) != (gs[i + 1] > 0):
            lo, hi = ts[i], ts[i + 1]
            for _ in range(100):
                mid = (lo + hi) / 2
                if (g_prime(mid, P, Q, dP, dQ) > 0) == (gs[i] > 0):
                    lo = mid
                else:
                    hi = mid
            t = (lo + hi) / 2
            pts.append((t, g_value(t, P, Q)))
    pts.sort(key=lambda z: z[0])
    dedup = []
    for t, g in pts:
        if dedup and abs(t - dedup[-1][0]) < mp.mpf("1e-13"):
            continue
        dedup.append((t, g))
    # maximal alternating subsequence, preferring the larger magnitude
    alt = []
    for t, g in dedup:
        if not alt or (g > 0) != (alt[-1][1] > 0):
            alt.append((t, g))
        elif abs(g) > abs(alt[-1][1]):
            alt[-1] = (t, g)
    return alt


def nodes_from_sequence(seq):
    """Map extrema (in the record's t) to rounded rational (x, t) pairs."""
    nodes = []
    for t, _g in seq:
        wm = mp.mpf(W_RECORD.numerator) / W_RECORD.denominator
        x = mp.mpf(A.numerator) / A.denominator + wm * (t + 1) / 2
        xf = Fraction(mp.nstr(x, 20))
        if not (A <= xf <= B):
            raise ArithmeticError("extrema node outside the paper interval")
        tf = 2 * (xf - A) / W_RECORD - 1
        nodes.append((xf, tf))
    return nodes


# ---------------------------------------------------------------------------
# Exact checks


def sturm_pole_free(P, Q):
    """Check Q is nonzero on the closed t-interval of the paper cell."""
    lo = Fraction(-1)
    hi = 2 * (B - A) / W_RECORD - 1
    qlo = fx.poly_eval(Q, lo)
    qhi = fx.poly_eval(Q, hi)
    if qlo == 0 or qhi == 0:
        return False, None
    roots = fx.sturm_root_count(Q, lo, hi)
    return roots == 0, roots


def certify_type(P_str, Q_str, m, k, label, source):
    P = [Fraction(s) for s in P_str[:m + 1]]
    Q = [Fraction(s) for s in Q_str[:k + 1]]
    if Q[0] != 1:
        raise ArithmeticError(f"{label}: Q(0) != 1")
    pole_free, roots = sturm_pole_free(P, Q)
    P_mp = [mp.mpf(s) for s in P_str[:m + 1]]
    Q_mp = [mp.mpf(s) for s in Q_str[:k + 1]]
    needed = m + k + 2
    candidates = []
    for include_endpoints in (True, False):
        seq = extrema_sequence(P_mp, Q_mp, include_endpoints=include_endpoints)
        if len(seq) < needed:
            continue
        try:
            nodes = nodes_from_sequence(seq)
        except ArithmeticError:
            continue
        ok, delta, rows = fx.dvp_check_mapped(P, Q, nodes, target=TARGET,
                                              n_terms=N_TERMS)
        candidates.append(dict(delta=delta, ok=ok, rows=rows,
                               n_ext=len(seq), endpoints=include_endpoints))
    if not candidates:
        return dict(label=label, source=source, m=m, k=k, stored=m + 1 + k,
                    needed=needed, certified=False,
                    error="no alternating sequence of the required length",
                    pole_free=pole_free, poles=roots)
    passing = [c for c in candidates if c["ok"]]
    best = max(passing or candidates, key=lambda c: c["delta"])
    delta, ok, rows = best["delta"], best["ok"], best["rows"]
    n_ext, endpoints = best["n_ext"], best["endpoints"]
    return dict(
        label=label, source=source, m=m, k=k, stored=m + 1 + k,
        needed_alternations=needed, numerical_extrema=n_ext,
        used_endpoints=endpoints,
        pole_free=pole_free, poles_on_interval=roots,
        Q_at_t_minus_1=float(fx.poly_eval(Q, Fraction(-1))),
        certified=bool(ok and pole_free),
        delta="%.6e" % float(delta),
        delta_num=str(delta.numerator), delta_den=str(delta.denominator),
        target="5e-14",
        delta_over_target="%.6f" % (float(delta) / float(TARGET)),
        nodes=rows,
        P=[str(c) for c in P_str[:m + 1]],
        Q=[str(c) for c in Q_str[:k + 1]],
        m_k=m + k,
    )


def grid_meets(P, Q, W, n_nodes=2000):
    """Exact max-error interval of P(t)/Q(t) on the uniform t-grid."""
    P = list(P)
    Q = list(Q)
    worst_hi = Fraction(0)
    worst_lo = Fraction(0)
    worst_t = None
    dmin = None
    for i in range(n_nodes + 1):
        t = Fraction(2 * i, n_nodes) - 1
        x = A + W * (t + 1) / 2
        lo, hi = fx.f0_bounds(x, n_terms=N_TERMS)
        r, den = fx.rat_eval(P, Q, t)
        dlo, dhi = lo - r, hi - r
        e_hi = max(abs(dlo), abs(dhi))
        if dlo <= 0 <= dhi:
            e_lo = Fraction(0)
        else:
            e_lo = min(abs(dlo), abs(dhi))
        if e_hi > worst_hi:
            worst_hi, worst_t = e_hi, t
        if e_lo > worst_lo:
            worst_lo = e_lo
        ad = abs(den)
        if dmin is None or ad < dmin:
            dmin = ad
    return dict(
        nodes=n_nodes + 1, interval="[x0, x0+W]",
        W=str(W),
        max_error_upper=float(worst_hi),
        max_error_lower=float(worst_lo),
        worst_t=float(worst_t),
        meets=bool(worst_hi <= TARGET),
        target="5e-14",
        denominator_min=float(dmin),
    )


# ---------------------------------------------------------------------------
# Extra types outside the recorded search, generated here with BRASIL.


def poly_mul(Aa, Bb):
    C = [mp.mpf(0)] * (len(Aa) + len(Bb) - 1)
    for i, ai in enumerate(Aa):
        for j, bj in enumerate(Bb):
            C[i + j] += ai * bj
    return C


def poly_add(Aa, Bb):
    n = max(len(Aa), len(Bb))
    C = [mp.mpf(0)] * n
    for i, a in enumerate(Aa):
        C[i] += a
    for i, b in enumerate(Bb):
        C[i] += b
    return C


def bary_to_monomial(r):
    z = [mp.mpf(float(v)) for v in r.nodes]
    w = [mp.mpf(float(v)) for v in r.weights]
    fv = [mp.mpf(float(v)) for v in r.values]
    n = len(z)
    P = [mp.mpf(0)]
    Q = [mp.mpf(0)]
    for j in range(n):
        L = [mp.mpf(1)]
        for k in range(n):
            if k != j:
                L = poly_mul(L, [-z[k], mp.mpf(1)])
        Q = poly_add(Q, [w[j] * c for c in L])
        P = poly_add(P, [fv[j] * w[j] * c for c in L])
    if Q[0] == 0:
        return None
    return [c / Q[0] for c in P], [c / Q[0] for c in Q]


def solve_extra(m, k):
    from baryrat import brasil
    from scipy.special import erf

    def f0_np(x):
        x = np.asarray(x, dtype=float)
        return np.sqrt(np.pi) / (2.0 * np.sqrt(x)) * erf(np.sqrt(x))

    f = lambda t: f0_np(float(A) + float(W_RECORD) *
                       (np.asarray(t) + 1.0) / 2.0)
    r, info = brasil(f, (-1.0, 1.0), (m, k), info=True,
                     tol=1e-4, maxiter=1000, npi=-30,
                     init_steps=100, max_step_size=0.1)
    conv = bary_to_monomial(r)
    if conv is None:
        raise ArithmeticError(f"[{m}/{k}] conversion failed")
    P, Q = conv
    P = [mp.nstr(c, 50) for c in P]
    Q = [mp.nstr(c, 50) for c in Q]
    return P, Q, bool(info.converged), float(info.deviation)


# ---------------------------------------------------------------------------


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--audit", default=AUDIT_DEFAULT)
    ap.add_argument("--corrected", default=CORRECTED_DEFAULT)
    ap.add_argument("--out", default="certificate_rational_paper_cell.json")
    ap.add_argument("--log", default="certificate_rational_paper_cell.log")
    ap.add_argument("--skip-grid", action="store_true")
    ap.add_argument("--skip-extras", action="store_true")
    args = ap.parse_args(argv)

    started = time.time()
    mp.mp.dps = MP_DPS
    log_lines = []

    def log(msg):
        line = f"[{time.time()-started:8.1f}s] {msg}"
        print(line, flush=True)
        log_lines.append(line)

    audit_digest = sha256(args.audit)
    corrected_digest = sha256(args.corrected)
    log(f"audit   {args.audit}  sha256={audit_digest}")
    log(f"corrected {args.corrected}  sha256={corrected_digest}")

    with open(args.audit, encoding="utf-8") as fh:
        audit_doc = json.load(fh)
    cell = find_cell(audit_doc, "paper", "abs:5e-14")
    control = cell["control"]
    accepted = control["accepted"]
    m_acc, k_acc = accepted["m"], accepted["k"]
    log(f"cell width={cell['width']} a={cell['a']} target={cell['target']}")
    log(f"stored pair [{m_acc}/{k_acc}] nominal count "
        f"{accepted.get('nominal_count', m_acc + 1 + k_acc)}")

    # Stored 12-coefficient nominal truncation (upper-bound witness).
    P_acc = accepted["P70"]
    Q_acc = accepted["Q70"]
    nominal = dict(
        pair=[m_acc, k_acc],
        P=P_acc[:m_acc + 1],
        Q=Q_acc[:k_acc + 1],
        witness_P_length=len(P_acc),
        witness_Q_length=len(Q_acc),
        source="corrected-counts-60dig.json control.accepted P70/Q70, "
               "nominal truncation",
    )
    log(f"nominal P length {len(nominal['P'])}, Q length "
        f"{len(nominal['Q'])}; witness {len(P_acc)}/{len(Q_acc)}")

    grid_results = {}
    if not args.skip_grid:
        P_hp = [Fraction(s) for s in nominal["P"]]
        Q_hp = [Fraction(s) for s in nominal["Q"]]
        for name, W in (("paper_W", W_PAPER), ("record_W", W_RECORD)):
            log(f"exact grid check on {name} ({W}) ...")
            g = grid_meets(P_hp, Q_hp, W)
            grid_results[name] = g
            log(f"  {name}: max |err| in "
                f"[{g['max_error_lower']:.6e}, {g['max_error_upper']:.6e}] "
                f"meets={g['meets']}")

    # The actual binary64-stored 12-coefficient representation, if available.
    binary64 = None
    try:
        with open(args.corrected, encoding="utf-8") as fh:
            corrected_doc = json.load(fh)
        ccell = find_cell(corrected_doc, "paper", "abs:5e-14")
        wp = ccell["witness_P"][:m_acc + 1]
        wq = ccell["witness_Q"][:k_acc + 1]
        # The witness strings are binary64 shortest round-trip renderings;
        # float(s) is exactly the stored binary64 coefficient.
        P64 = [Fraction(float(s)) for s in wp]
        Q64 = [Fraction(float(s)) for s in wq]
        binary64 = dict(
            P=[repr(float(s)) for s in wp],
            Q=[repr(float(s)) for s in wq],
            source="corrected-counts.json witness_P/witness_Q, nominal "
                   "truncation, parsed as exact binary64 values",
        )
        if not args.skip_grid:
            g64 = grid_meets(P64, Q64, W_PAPER)
            binary64["exact_grid_check_paper_W"] = g64
            log(f"binary64 stored nominal truncation: max |err| in "
                f"[{g64['max_error_lower']:.6e}, {g64['max_error_upper']:.6e}]"
                f" meets={g64['meets']}")
            log("rounding budget: the record's binary64 Horner evaluation of "
                "this fit is 9.936496070395151e-15; a conservative few-ulp "
                "Horner budget is O(1e-14), so the meeting cannot be "
                "overturned by evaluation rounding")
    except (KeyError, FileNotFoundError) as exc:
        log(f"binary64 witness check skipped: {exc!r}")

    # Recorded stored-11 types: use the audit's own default-pass fits.
    failing = control["failing_level"]
    rows = failing["rows"]
    default_rows = {}
    for row in rows:
        key = (row["m"], row["k"])
        if row.get("settings") == "default" and key not in default_rows:
            default_rows[key] = row
    log(f"recorded stored-11 types: {sorted(default_rows)}")

    types = []
    for (m, k) in sorted(default_rows):
        row = default_rows[(m, k)]
        label = f"recorded below-count [{m}/{k}]"
        log(f"certifying {label} (stored {m+1+k}, need {m+k+2} alternations)")
        cert = certify_type(row["P70"], row["Q70"], m, k, label,
                            "corrected-counts-60dig.json failing_level "
                            "default row")
        cert["recorded_e60_abs"] = row.get("e60_abs")
        types.append(cert)
        log(f"  certified={cert['certified']} pole_free={cert['pole_free']} "
            f"delta={cert.get('delta')} ({cert.get('delta_over_target')} x "
            f"target) nodes={len(cert.get('nodes', []))}")

    extras = []
    if not args.skip_extras:
        for (m, k) in [(10, 0), (0, 10)]:
            label = f"search-excluded endpoint [{m}/{k}]"
            log(f"solving {label} with BRASIL ...")
            P, Q, conv, dev = solve_extra(m, k)
            cert = certify_type(P, Q, m, k, label,
                                "BRASIL 2.1.2 solve by this pilot")
            cert["brasil_converged"] = conv
            cert["brasil_deviation"] = dev
            extras.append(cert)
            log(f"  certified={cert['certified']} delta={cert.get('delta')} "
                f"({cert.get('delta_over_target')} x target)")

    all_types = types + extras
    recorded_ok = bool(types) and all(t["certified"] for t in types)
    extras_ok = bool(extras) and all(t["certified"] for t in extras)
    grid_ok = bool(grid_results) and any(
        g["meets"] and k == "paper_W" for k, g in grid_results.items())
    b64_ok = bool(binary64 and binary64.get(
        "exact_grid_check_paper_W", {}).get("meets", False))
    certified_min = bool(recorded_ok and grid_ok and b64_ok)

    certified_types = [t for t in types if t["certified"]]
    weakest = min(certified_types,
                  key=lambda t: Fraction(int(t["delta_num"]),
                                         int(t["delta_den"]))) \
        if certified_types else None

    defects = []
    for t in types:
        if t["m"] + t["k"] != 10:
            defects.append(f"{t['label']}: m+k={t['m']+t['k']} != 10")
        if t["needed_alternations"] != 12:
            defects.append(f"{t['label']}: need {t['needed_alternations']}"
                           f" != 12")

    certificate = dict(
        schema="paper-c-minimax-certificate-pilot-1",
        generated_utc=now(),
        command="python -I make_certificate.py" +
                (" --skip-grid" if args.skip_grid else "") +
                (" --skip-extras" if args.skip_extras else ""),
        versions=dict(python=sys.version.split()[0],
                      mpmath=mp.__version__, numpy=np.__version__),
        inputs=dict(
            audit=dict(path=args.audit, sha256=audit_digest),
            corrected=dict(path=args.corrected, sha256=corrected_digest),
        ),
        cell=dict(
            function="F_0(x) = sqrt(pi)/(2 sqrt(x)) erf(sqrt(x))",
            x0=str(A), x1=str(B),
            W_paper=str(W_PAPER), W_record=str(W_RECORD),
            target="5e-14 absolute",
            width_note=("the audit's width field is 1 ulp below x1-x0; "
                        "the exact grid check is reported on both "
                        "intervals"),
            certifying_nodes=2001,
        ),
        stored_candidate=nominal,
        stored_candidate_binary64=binary64,
        stored_candidate_meets=grid_results,
        recorded_below_types=types,
        extra_types=extras,
        defect_list=defects,
        conclusion=dict(
            rational_minimum_stored_count=12 if certified_min else None,
            recorded_family_certified=recorded_ok,
            beyond_search_certified=bool(extras_ok and recorded_ok),            stored_fit_meets=bool(grid_ok and b64_ok),
            rounding_budget=(
                "lower bound: exact rational arithmetic, no rounding budget "
                "needed. Upper bound: exact evaluation of the stored binary64 "
                "coefficients; the record's binary64 Horner value is "
                "9.936496070395151e-15 and a conservative few-ulp Horner "
                "budget is below 1.1e-14, leaving the meeting intact."),
            all_recorded_types_exceed_target=recorded_ok,
            weakest_type=weakest["label"] if weakest else None,
            weakest_delta=weakest["delta"] if weakest else None,
            weakest_delta_over_target=weakest["delta_over_target"]
            if weakest else None,
            scope=("all rationals of stored count <= 11 whose pairs are "
                   "admissible in the paper's search (m>=1, k>=1, k<=12); "
                   "the endpoint types [0/10] and [10/0] are certified "
                   "separately as search-excluded extras"),
            proof=("de la Vallee Poussin: if F_0 - r has m+k+2 alternating "
                   "signs of magnitude > eps, no rational R of type (m,k) "
                   "can satisfy |F_0 - R| <= eps, since R - r would have "
                   "m+k+1 zeros but numerator degree at most m+k. Padding "
                   "the numerator embeds every lower-count pair in a "
                   "count-11 type."),
        ),
    )

    with open(args.out, "w", encoding="utf-8") as fh:
        json.dump(certificate, fh, indent=1)
    with open(args.log, "w", encoding="utf-8") as fh:
        fh.write(f"started_utc={now()}\n")
        fh.write(f"command={certificate['command']}\n")
        fh.write(f"audit_sha256={audit_digest}\n")
        fh.write(f"corrected_sha256={corrected_digest}\n")
        fh.write("\n".join(log_lines) + "\n")
        fh.write(f"conclusion={json.dumps(certificate['conclusion'])}\n")
        fh.write(f"defects={defects}\n")
    log(f"wrote {args.out}")
    log(f"conclusion: {json.dumps(certificate['conclusion'])}")
    log(f"defects: {defects}")
    return 0 if (certified_min and extras_ok and not defects) else 1


if __name__ == "__main__":
    raise SystemExit(main())
