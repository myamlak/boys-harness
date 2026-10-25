#!/usr/bin/env python3
"""Explore candidate monotone freezes and find the one that reproduces the
manuscript's 21 commits / 21 right / 46 lost.  Read-only."""

import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import surrogate_control as SC  # noqa: E402
from mpmath import mp, mpf  # noqa: E402

mp.dps = 70
A = SC.A
NMAX = SC.NMAX


def run(label, fun, entries, winner, margin=mpf(1)):
    rec_abs = SC.predicted(entries, fun, "absolute")
    rec_rel = SC.predicted(entries, fun, "relative")
    table = SC.all_rows(*rec_abs, *rec_rel, winner, margin)
    total, per = SC.score_all(table)
    print("%-52s commits %3d right %3d wrong %3d   " % (label, *total), end="")
    print(
        " | ".join("%s/%s %d/%d/%d" % (k[0], k[1][:3], *per[k]) for k in sorted(per))
    )
    return table


def main():
    reference, cells, error0, winner, *_ = SC.build_measured()
    entries, population, model_from = SC.build_cells(reference, cells, error0)

    def nmin_ceil(x):
        return min(NMAX, max(0, math.ceil(x - 0.5)))

    def nmin_floor(x):
        return min(NMAX, max(0, math.floor(x + 0.5)))

    def nmin_round(x):
        return min(NMAX, max(0, round(x - 0.5)))

    def freeze_after(nmin_of):
        def fun(n, x):
            m = nmin_of(x)
            if n <= m:
                return A.amplification(n, x)
            return A.amplification(m, x)

        return fun

    def freeze_all(nmin_of):
        def fun(n, x):
            return A.amplification(nmin_of(x), x)

        return fun

    def freeze_all_n0(nmin_of):
        def fun(n, x):
            if n == 0:
                return mpf(1)
            return A.amplification(nmin_of(x), x)

        return fun

    def unit(n, x):
        return mpf(1)

    print("population A1 %d A2 %d" % (population["A1"], population["A2"]))
    print("margin 1:")
    run("V0 true model", A.amplification, entries, winner)
    run("V1 freeze after nmin (ceil x-0.5), falling kept", freeze_after(nmin_ceil), entries, winner)
    run("V2 freeze at min for every order (ceil x-0.5)", freeze_all(nmin_ceil), entries, winner)
    run("V2n0 freeze at min for n>=1, A0=1", freeze_all_n0(nmin_ceil), entries, winner)
    run("V2f freeze at min every order (floor x+0.5)", freeze_all(nmin_floor), entries, winner)
    run("V2r freeze at min every order (round x-0.5)", freeze_all(nmin_round), entries, winner)
    run("V1f freeze after (floor x+0.5)", freeze_after(nmin_floor), entries, winner)
    run("V3 A_n = 1", unit, entries, winner)
    return 0


if __name__ == "__main__":
    sys.exit(main())
