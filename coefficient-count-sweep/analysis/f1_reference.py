#!/usr/bin/env python3
"""The F_1 reference for the second instance of the tail class.

Paper C's counts are measured on F_0.  The tail-class proposition of
main.tex (section "The tail-class proposition") admits any algebraic tail
f = c x^(-alpha) g(x) with g analytic and non-vanishing on the band, and it
states the Chebyshev rate log(1/kappa) as exponent-independent.  F_1 is a
second instance at alpha = 3/2:

    F_1(x) = (F_0(x) - exp(-x)) / (2 x)          (the n = 0 -> 1 recurrence)
           = sum_{m>=0} (-1)^m x^m / (m! (2m+3)) (entire, from the series)
           ~ (sqrt(pi)/4) x^(-3/2) g(x),  g -> 1     (the algebraic tail)

The recurrence is exact for the Boys functions,
F_{n+1}(x) = ((2n+1) F_n(x) - exp(-x)) / (2x), and F_1(x) is also the exact
integral int_0^1 t^2 exp(-x t^2) dt.  This module supplies

  * f1_reference(g, xs)   the 60-digit reference on a grid, from the paper's
                          own Boys routine boys_ref(0, x, RAT_TAIL_FLOOR,
                          RAT_TERMS) followed by the recurrence;
  * f1_cheb_fit(g, a, b, deg)
                          the Chebyshev interpolation of F_1 at the same
                          nodes the generator's cheb_coeffs uses, through the
                          generator's own cheb_coeffs_of;
  * a validation CLI that checks the reference against the exact F_0
                          identity, against the high-precision integral, and
                          against the entire series.

Nothing here reimplements a solver, an evaluator or an acceptance test: the
driver monkey-patches count_sweep.cheb_fit onto f1_cheb_fit and hands the
F_1 reference to count_sweep.cheb_count and count_sweep.rational_count, which
are the paper's own search, acceptance and evaluation code.
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import mpmath as mp

PAPER = Path(__file__).resolve().parent
sys.path.insert(0, str(PAPER))

import count_sweep as cs  # noqa: E402

X0_TEXT = "11.899848152108484"
# The sample grid for the validation: inside both bands, at both left edges,
# at the widest grid point of the second edge, and two small-x stress points.
SAMPLE_X = [
    "0.125", "0.25", "0.5", "1", "2",
    X0_TEXT,
    "23.799696304216968",
    "43.899848152108487",
    "87.799696304216968",
    "200",
]


def f1_from_f0(g, x, tail_floor=None, terms=None):
    """F_1(x) from the recurrence, F_0 from the paper's Boys routine.

    ``tail_floor``/``terms`` default to the fit path's precision
    (RAT_TAIL_FLOOR=1e-55, RAT_TERMS=900), which is the precision
    count_sweep.reference certifies against.  The subtraction loses digits
    only by the ratio exp(-x)/F_0(x), 2.7e-5 at x0 and far smaller outward, so
    the reference keeps more than fifty significant digits on the band.
    """
    floor = g.RAT_TAIL_FLOOR if tail_floor is None else tail_floor
    nterms = g.RAT_TERMS if terms is None else terms
    f0 = g.boys_ref(0, x, floor, nterms)
    return (f0 - g.mp.exp(-x)) / (2 * x)


def f1_from_f0_fit(g, x):
    """F_1 at the generator's *fit-path* Boys precision (REF_TAIL_FLOOR).

    cheb_coeffs(f) computes its interpolation values with boys_ref's default
    arguments (REF_TAIL_FLOOR=1e-26, REF_TERMS=300), not the 1e-55 reference
    the acceptance test uses.  The Chebyshev coefficients are therefore
    formed from this value so the fit path is the one the F_0 sweep used;
    the acceptance test still sees the 1e-55 reference.  Twenty-six digits
    against a 1e-14 target is not a boundary.
    """
    return f1_from_f0(g, x, g.REF_TAIL_FLOOR, g.REF_TERMS)


def f1_reference(g, xs):
    """The 60-digit F_1 reference on a grid, exactly as count_sweep.reference."""
    with g.mp.workdps(g.RAT_DPS):
        return [f1_from_f0(g, x) for x in xs]


def f1_cheb_fit(g, a, b, deg):
    """The Chebyshev coefficients of F_1 at the generator's interpolation nodes.

    Drop-in replacement for count_sweep.cheb_fit (which calls
    g.cheb_coeffs(0, a, b, deg) and is therefore tied to F_0).  The nodes,
    the cosine sums and the rounding are the generator's cheb_coeffs_of, the
    function the generator already carries for "fits whose target is not F
    itself".
    """
    with g.mp.workdps(g.RAT_DPS):
        cs_coeffs = g.cheb_coeffs_of(lambda x: f1_from_f0_fit(g, x), a, b, deg)
    return [float(c) for c in cs_coeffs]


# ---------------------------------------------------------------------------
# validation


def f0_identity(g, x):
    """F_0 from the exact closed form, (sqrt(pi)/(2 sqrt(x))) erf(sqrt(x))."""
    return g.mp.sqrt(g.mp.pi) * g.mp.erf(g.mp.sqrt(x)) / (2 * g.mp.sqrt(x))


def f1_integral(g, x):
    """F_1 from its definition, int_0^1 t^2 exp(-x t^2) dt."""
    return g.mp.quad(lambda t: t ** 2 * g.mp.exp(-x * t ** 2), [0, 1])


def f1_series(g, x):
    """F_1 from its entire Taylor series, sum (-1)^m x^m/(m!(2m+3))."""
    x = g.mp.mpf(x)
    total = g.mp.mpf(0)
    term = g.mp.mpf(1) / 3
    m = 0
    scale = abs(term)
    while True:
        total += term
        m += 1
        term = -term * x * (2 * m + 1) / (m * (2 * m + 3))
        if abs(term) < scale * g.mp.mpf("1e-70") and m > 2 * int(x) + 20:
            break
        if m > 10000:
            raise RuntimeError(f"series did not converge at x={g.mp.nstr(x, 12)}")
    return total


def recurrence_check(g, x):
    """The staircase: ((2n+1)F_n - exp(-x))/(2x) against boys_ref(n+1, x)."""
    out = {}
    with g.mp.workdps(g.RAT_DPS):
        f = g.boys_ref(0, x, g.RAT_TAIL_FLOOR, g.RAT_TERMS)
        ex = g.mp.exp(-x)
        for n in (0, 1, 2):
            step = ((2 * n + 1) * f - ex) / (2 * x)
            direct = g.boys_ref(n + 1, x, g.RAT_TAIL_FLOOR, g.RAT_TERMS)
            out[str(n)] = {
                "relative_difference": g.mp.nstr(abs(step - direct) / direct, 10),
            }
            f = direct
    return out


def digits_agree(value, reference):
    """The number of leading decimal digits at which two mpf values agree."""
    if value == reference:
        return 99
    rel = abs(value - reference) / max(abs(reference), mp.mpf("1e-300"))
    if rel == 0:
        return 99
    return float(-mp.log10(rel))


def validate(boys: Path, out_path: Path | None = None) -> dict:
    g = cs.load_generator(boys)
    record = {
        "kernel_root": str(boys),
        "generator_sha256": cs.sha256(boys / cs.GENERATOR_REL),
        "header_sha256": cs.sha256(boys / cs.HEADER_REL),
        "RAT_DPS": g.RAT_DPS,
        "RAT_TAIL_FLOOR": g.mp.nstr(g.RAT_TAIL_FLOOR, 3),
        "RAT_TERMS": g.RAT_TERMS,
        "identity": (
            "F_1(x) = (F_0(x) - exp(-x))/(2x); F_0 from "
            "boys_ref(0, x, RAT_TAIL_FLOOR, RAT_TERMS) at RAT_DPS; the exact "
            "closed form and the integral are separate witnesses"),
        "samples": {},
    }
    for text in SAMPLE_X:
        with g.mp.workdps(g.RAT_DPS):
            x = g.mpf(text)
            f0_ref = g.boys_ref(0, x, g.RAT_TAIL_FLOOR, g.RAT_TERMS)
            f0_exact = f0_identity(g, x)
            f1_ref = f1_from_f0(g, x)
            with g.mp.workdps(g.RAT_DPS + 20):
                # The witnesses run at 80 digits so their own rounding does
                # not enter the comparison at 60.
                xx = g.mpf(text)
                integral = f1_integral(g, xx)
            # The alternating Taylor series has peak term ~exp(x)/sqrt(x), so
            # it cancels about x*log10(e) digits; run it where that is free of
            # the working precision.  This is a property of the witness, not
            # of the reference, which the integral certifies everywhere.
            series_dps = g.RAT_DPS + 20 + int(0.5 * float(x))
            with g.mp.workdps(series_dps):
                series = f1_series(g, g.mpf(text))
            tail = (x ** g.mpf("1.5")) * f1_ref / (g.mp.sqrt(g.mp.pi) / 4)
            entry = {
                "x": text,
                "F0_boys_ref": g.mp.nstr(f0_ref, 25),
                "F0_exact": g.mp.nstr(f0_exact, 25),
                "F0_boys_ref_vs_exact_digits": digits_agree(f0_ref, f0_exact),
                "F1_reference": g.mp.nstr(f1_ref, 25),
                "F1_integral": g.mp.nstr(integral, 25),
                "F1_reference_vs_integral_digits": digits_agree(f1_ref, integral),
                "F1_series": g.mp.nstr(series, 25),
                "F1_series_dps": series_dps,
                "F1_reference_vs_series_digits": digits_agree(f1_ref, series),
                "x_pow_3_2_F1_over_sqrtpi_over_4":
                    g.mp.nstr(tail, 20),
                "recurrence_vs_boys_ref_digits": recurrence_check(g, x),
            }
            record["samples"][text] = entry
            print(f"x={text:>22}  F0 ref~exact {entry['F0_boys_ref_vs_exact_digits']:5.1f} digits"
                  f"  F1 ref~integral {entry['F1_reference_vs_integral_digits']:5.1f} digits"
                  f"  F1 ref~series {entry['F1_reference_vs_series_digits']:5.1f} digits"
                  f"  x^(3/2)F1/(sqrt(pi)/4) = {entry['x_pow_3_2_F1_over_sqrtpi_over_4']}",
                  flush=True)
    worst_integral = min(e["F1_reference_vs_integral_digits"]
                         for e in record["samples"].values())
    worst_series = min(e["F1_reference_vs_series_digits"]
                       for e in record["samples"].values())
    worst_f0 = min(e["F0_boys_ref_vs_exact_digits"]
                    for e in record["samples"].values())
    record["worst_digits"] = {
        "F0_boys_ref_vs_exact": worst_f0,
        "F1_reference_vs_integral": worst_integral,
        "F1_reference_vs_series": worst_series,
    }
    record["verdict"] = (
        "the F_1 reference agrees with the high-precision integral and with "
        f"the entire series to at least {min(worst_integral, worst_series):.0f} "
        f"digits, and the F_0 reference it is built on agrees with the exact "
        f"closed form to at least {worst_f0:.0f} digits")
    print("\n" + record["verdict"], flush=True)
    if out_path is not None:
        out_path.parent.mkdir(parents=True, exist_ok=True)
        out_path.write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8")
        print(f"wrote {out_path}", flush=True)
    return record


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--boys-root", default=r"C:\Users\myaml\dev\_mainonly")
    parser.add_argument("--out", default="")
    args = parser.parse_args(argv)
    started = time.perf_counter()
    out = Path(args.out) if args.out else None
    validate(Path(args.boys_root), out)
    print(f"validation wall time: {time.perf_counter() - started:.1f} s", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
