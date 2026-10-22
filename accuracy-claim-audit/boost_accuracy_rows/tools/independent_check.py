"""Ask an outside implementation whether the gate's reference values are what it says.

The gate measures against boost::multiprecision::cpp_bin_float_100 evaluated by the same library
it is measuring.  A reference that shares an implementation with the thing under test can be
wrong in the same way, so this recomputes the same points with mpmath -- a different codebase,
by a different author -- at 60 significant digits, and reports the disagreement in units of
double epsilon.

It also anchors two points that are exact by definition, and one against a constant that comes
from outside either library:

    tgamma(5) = 4! = 24 exactly
    lgamma(2) = log(1!) = 0 exactly
    digamma(1) = -euler_gamma, with the constant taken from mpmath

Usage: python tools/independent_check.py claims/reference_dump.csv
"""

import csv
import sys
from pathlib import Path

import mpmath as mp

mp.mp.dps = 60

EVAL = {
    "erf": lambda nu, x: mp.erf(x),
    "erfc": lambda nu, x: mp.erfc(x),
    "tgamma": lambda nu, x: mp.gamma(x),
    "lgamma": lambda nu, x: mp.loggamma(x),
    "digamma": lambda nu, x: mp.digamma(x),
    "cyl_bessel_j": lambda nu, x: mp.besselj(nu, x),
    "cyl_bessel_i": lambda nu, x: mp.besseli(nu, x),
    "expint": lambda nu, x: mp.ei(x),
    "airy_ai": lambda nu, x: mp.airyai(x),
    "airy_bi": lambda nu, x: mp.airybi(x),
}

EPS = mp.mpf(2) ** -52


def relative_difference(a, b):
    """The library's own metric: max of the two one-sided relative differences."""
    if a == b:
        return mp.mpf(0)
    return max(abs((a - b) / a), abs((a - b) / b))


def main() -> int:
    dump = list(csv.DictReader(Path(sys.argv[1]).open(encoding="utf-8")))
    by_function = {}
    for row in dump:
        by_function.setdefault(row["function"], []).append(row)

    worst_agreement = mp.mpf(0)
    worst_at = None
    for function, rows in sorted(by_function.items()):
        step = max(1, len(rows) // 8)
        picks = rows[::step]
        if rows[-1] not in picks:
            picks.append(rows[-1])
        worst = mp.mpf(0)
        worst_row = None
        for row in picks:
            nu, x = mp.mpf(row["order"]), mp.mpf(row["x"])
            mine = mp.mpf(row["reference"])
            theirs = EVAL[function](nu, x)
            diff = relative_difference(mine, theirs) / EPS
            if diff > worst:
                worst, worst_row = diff, (nu, x, mine, theirs)
        nu, x, mine, theirs = worst_row
        print(f"{function:14s} {len(picks):3d} points, worst gate-vs-mpmath "
              f"{mp.nstr(worst, 6):>10s} eps at x={mp.nstr(x, 12)}")
        if worst > worst_agreement:
            worst_agreement, worst_at = worst, (function, x)
    print()
    print(f"worst disagreement anywhere: {mp.nstr(worst_agreement, 6)} eps "
          f"({worst_at[0]} at x={mp.nstr(worst_at[1], 12)})")
    print(f"double epsilon is 1 eps, so the reference is independent to "
          f"{mp.nstr(1 / worst_agreement, 3)}x the resolution of the comparison"
          if worst_agreement else "")
    print()

    print("exact anchors")
    t5 = mp.gamma(5)
    l2 = mp.loggamma(2)
    d1 = mp.digamma(1)
    print(f"  mpmath tgamma(5)  = {mp.nstr(t5, 20)}   exact 24: {t5 == 24}")
    print(f"  mpmath lgamma(2)  = {mp.nstr(l2, 20)}   exact 0 : {l2 == 0}")
    print(f"  mpmath digamma(1) = {mp.nstr(d1, 20)}   -euler_gamma: "
          f"{mp.nstr(d1 + mp.euler, 3)} off")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
