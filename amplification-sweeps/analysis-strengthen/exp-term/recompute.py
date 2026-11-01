#!/usr/bin/env python3
"""Every ablation and population count that depends on the model's committed set.

The surrogate control is imported from the paper-tree copy beside this file; it
imports the regenerated analyze.py, so its arithmetic term carries the
exponential contribution.  This driver does not modify it: it calls its
builders and scores the results itself, at the regenerated model's own
committed count rather than the previous 67.

Run without python -I (user-site mpmath).
"""

import csv
import math
import sys
from pathlib import Path

from mpmath import mp, mpf

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE / "surrogate"))
import surrogate_control as S  # noqa: E402

mp.dps = 70
NMAX = S.NMAX
BANDS = ("A1", "A2")
METRICS = ("absolute", "relative")
THETA0 = 1.0855252345349333
RESULTS = S.TREE / "results"


def table_set(table):
    return {
        (key[0], key[1], row["n"])
        for key, rows in table.items()
        for row in rows
        if row["state"] != "abstain"
    }


def score(table):
    total, per = S.score_all(table)
    return total, per


def frontier_zero(steps):
    zero = sorted({s[2] for s in steps if s[4] == 0})
    return (zero[0] if zero else None, zero[-1] if zero else None, zero)


def step_for(steps, target):
    for i, s in enumerate(steps):
        if s[2] == target:
            return i, s
    return None, None


def main():
    reference, cells, error0, winner, *_ = S.build_measured()
    entries, population, model_from = S.build_cells(reference, cells, error0)
    print("population", population, "model_from", model_from)

    true_abs = S.predicted(entries, S.A.amplification, "absolute")
    true_rel = S.predicted(entries, S.A.amplification, "relative")
    true1 = S.all_rows(*true_abs, *true_rel, winner, mpf(1))
    (tc, tr, tw), per = score(true1)
    true_set = table_set(true1)
    print("TRUE margin 1: commits %d right %d wrong %d" % (tc, tr, tw))
    for key in sorted(per):
        print("  %s %-8s commits %d right %d wrong %d" % (key[0], key[1], *per[key]))
    true_steps, tac, taa, _ = S.frontier(*true_abs, *true_rel, winner)
    tz0, tz1, tzero = frontier_zero(true_steps)
    print("TRUE frontier zero-error commits: %d .. %d (%d distinct counts)"
          % (tz0, tz1, len(tzero)))
    # the error-free ceiling and the margin interval that attains it
    top = [s for s in true_steps if s[2] == tz1 and s[4] == 0]
    print("TRUE max zero-error step(s):", [(float(s[0]), float(s[1]), s[2], s[3], s[4]) for s in top])

    # ---- monotone surrogate ---------------------------------------------
    cache = {}

    def nmin_of(x):
        return min(NMAX, max(0, math.ceil(x - 0.5)))

    def frozen_min(n, x):
        m = nmin_of(x)
        if n <= m:
            return S.A.amplification(n, x)
        key = (m, x)
        if key not in cache:
            cache[key] = S.A.amplification(m, x)
        return cache[key]

    sur_abs = S.predicted(entries, frozen_min, "absolute")
    sur_rel = S.predicted(entries, frozen_min, "relative")
    sur1 = S.all_rows(*sur_abs, *sur_rel, winner, mpf(1))
    (sc, sr, sw), sper = score(sur1)
    sur_set = table_set(sur1)
    lost = true_set - sur_set
    gained = sur_set - true_set
    wrong_lost = [
        k for k in lost
        if next(r for r in true1[(k[0], k[1])] if r["n"] == k[2])["right"] is False
    ]
    wrong_gained = [
        k for k in gained
        if next(r for r in sur1[(k[0], k[1])] if r["n"] == k[2])["right"] is False
    ]
    print()
    print("MONOTONE SURROGATE margin 1: commits %d right %d wrong %d" % (sc, sr, sw))
    for key in sorted(sper):
        print("  %s %-8s commits %d right %d wrong %d" % (key[0], key[1], *sper[key]))
    print("  true commits %d, surrogate commits %d, lost %d, gained %d"
          % (len(true_set), len(sur_set), len(lost), len(gained)))
    print("  lost that were wrong for the true model: %d" % len(wrong_lost))
    print("  gained that are wrong for the surrogate: %d" % len(wrong_gained))
    sur_steps, sac, saa, _ = S.frontier(*sur_abs, *sur_rel, winner)
    sz0, sz1, szero = frontier_zero(sur_steps)
    print("  surrogate frontier zero-error commits: %d .. %d (%d distinct)"
          % (sz0, sz1, len(szero)))
    stop = [s for s in sur_steps if s[2] == sz1 and s[4] == 0]
    print("    max zero-error step(s):",
          [(float(s[0]), float(s[1]), s[2], s[3], s[4]) for s in stop])
    print("  shape worth (true ceiling - surrogate ceiling): %d" % (tz1 - sz1))
    for m in (mpf("0.94"), mpf(1), mpf("0.946")):
        tt, _ = score(S.all_rows(*true_abs, *true_rel, winner, m))
        ss, _ = score(S.all_rows(*sur_abs, *sur_rel, winner, m))
        print("  at margin %s: true %s, surrogate %s" % (mp.nstr(m, 5), tt, ss))
    target = tc
    idx, step = step_for(sur_steps, target)
    print("  matched target %d: %s" % (target, "attained" if step else "skipped"))
    if step:
        lo, hi, c, r, w = step
        m = (lo + hi) / 2
        t, _ = score(S.all_rows(*sur_abs, *sur_rel, winner, m))
        print("    margin interval [%s, %s), midpoint %s -> commits %d right %d wrong %d"
              % (mp.nstr(lo, 17), mp.nstr(hi, 17), mp.nstr(m, 17), *t))
        # neighbours
        for j in (idx - 1, idx + 1):
            if 0 <= j < len(sur_steps):
                lo2, hi2, c2, r2, w2 = sur_steps[j]
                print("    neighbour commits %d right %d wrong %d on [%s, %s)"
                      % (c2, r2, w2, mp.nstr(lo2, 17), mp.nstr(hi2, 17)))
        matched_margin = m
    else:
        matched_margin = mpf("0.946")

    # ---- frozen-first and unit controls ---------------------------------
    def frozen_first(n, x):
        return mpf(1) / (2 * x)

    def frozen_first_n0(n, x):
        return mpf(1) if n == 0 else mpf(1) / (2 * x)

    def unit(n, x):
        return mpf(1)

    for label, fun in (
        ("FROZEN-FIRST A_n=1/(2x), n=0 kept A_0=1", frozen_first_n0),
        ("FROZEN-FIRST A_n=1/(2x), n=0 also", frozen_first),
        ("UNIT A_n=1", unit),
    ):
        f_abs = S.predicted(entries, fun, "absolute")
        f_rel = S.predicted(entries, fun, "relative")
        f1 = S.all_rows(*f_abs, *f_rel, winner, mpf(1))
        (fc, fr, fw), fper = score(f1)
        fset = table_set(f1)
        flost = true_set - fset
        fgain = fset - true_set
        fwrong = [
            k for k in fgain
            if next(r for r in f1[(k[0], k[1])] if r["n"] == k[2])["right"] is False
        ]
        fsteps, fac, faa, _ = S.frontier(*f_abs, *f_rel, winner)
        fz0, fz1, fzero = frontier_zero(fsteps)
        print()
        print("%s margin 1: commits %d right %d wrong %d" % (label, fc, fr, fw))
        for key in sorted(fper):
            print("  %s %-8s commits %d right %d wrong %d" % (key[0], key[1], *fper[key]))
        print("  vs true set: lost %d, gained %d, gained-wrong %d"
              % (len(flost), len(fgain), len(fwrong)))
        if fz0 is None:
            best = min(fsteps, key=lambda s: (s[4], -s[2]))
            print("  zero-error commits: none; smallest wrong count %d at commits %d"
                  % (best[4], best[2]))
        else:
            print("  zero-error commits: %d .. %d (%d distinct)"
                  % (fz0, fz1, len(fzero)))
            top = [s for s in fsteps if s[2] == fz1 and s[4] == 0]
            print("    max zero-error step(s):",
                  [(float(s[0]), float(s[1]), s[2], s[3], s[4]) for s in top])

    # ---- population-restricted winners ----------------------------------
    winner_pop = S.winners_over(reference, cells, model_from)
    differing = [
        key for key in true_set
        if winner_pop[(key[0], key[1], key[2])] != winner[(key[0], key[1], key[2])]
    ]
    tpop = S.all_rows(*true_abs, *true_rel, winner_pop, mpf(1))
    (pc, pr, pw), _ = score(tpop)
    print()
    print("POPULATION winners x >= theta0, true model margin 1: commits %d right %d wrong %d"
          % (pc, pr, pw))
    print("  orders where the restricted winner differs among true commitments: %d" % len(differing))
    for key in sorted(differing):
        print("    %s %-8s n=%2d whole %-7s population %-7s"
              % (key[0], key[1], key[2],
                 winner[(key[0], key[1], key[2])],
                 winner_pop[(key[0], key[1], key[2])]))
    spop = S.all_rows(*sur_abs, *sur_rel, winner_pop, matched_margin)
    (spc, spr, spw), _ = score(spop)
    print("  surrogate at matched margin, population winners: commits %d right %d wrong %d"
          % (spc, spr, spw))

    # ---- measured winning cell below theta0 ------------------------------
    order = {}
    with open(RESULTS / "per_order.csv", newline="") as handle:
        for row in csv.DictReader(handle):
            order[(row["band"], int(row["n"]), row["method"])] = row
    below = []
    for band, metric, n in sorted(true_set):
        w = winner[(band, metric, n)]
        if w == "tie":
            continue
        at = "at_x_abs" if metric == "absolute" else "at_x_rel"
        x = float(order[(band, n, w)][at])
        if x < THETA0:
            below.append((band, metric, n, w, x))
    print()
    print("COMMITMENTS whose measured winning cell is below theta0: %d of %d"
          % (len(below), len(true_set)))
    for row in below:
        print("  %s %-8s n=%2d winner %-7s x=%.17g" % row)

    # first sweep, read from its prediction file
    first = list(csv.DictReader(open(HERE.parent / "out" / "first" / "prediction.csv")))
    fc = [r for r in first if r["state"] != "abstain"]
    fr = sum(1 for r in fc if r["right_where_committed"] == "1")
    print()
    print("FIRST SWEEP prediction: commits %d right %d wrong %d"
          % (len(fc), fr, len(fc) - fr))
    first_m0 = list(csv.DictReader(open(HERE.parent / "out" / "first_m0" / "prediction.csv")))
    f0 = [r for r in first_m0 if r["state"] != "abstain"]
    f0r = sum(1 for r in f0 if r["right_where_committed"] == "1")
    f0tie = sum(1 for r in f0 if r["measured_winner"] == "tie")
    print("FIRST SWEEP margin 0: commits %d right %d wrong %d ties %d"
          % (len(f0), f0r, len(f0) - f0r, f0tie))


if __name__ == "__main__":
    main()
