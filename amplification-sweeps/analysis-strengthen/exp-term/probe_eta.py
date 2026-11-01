#!/usr/bin/env python3
"""Compare three readings of 'one unit in the last place' for eta_E.

    const 2^-52   the largest relative ulp of a normal binary64
    const 2^-53   the unit roundoff of binary64
    ulp(value)    one ulp of the computed exp(-x) at each argument,
                  eta_E(x) = ulp(exp(-x))/exp(-x)

The model is score_exp.Model with eta_E passed in; eta_E = 0 reproduces the
recorded states (score-exp-run.log).  Prints the margin-1.0 score and the
rows whose state differs from the recorded one for each reading.
"""

import math
from mpmath import mp, mpf
import score_exp as S

mp.dps = 70


def eta_ulp(x):
    v = math.exp(-x)
    if v == 0.0:
        return mpf(2) ** -52
    return mpf(math.ulp(v)) / mpf(v)


def main():
    module = S.load_analysis()
    trees = {
        "second": (
            S.SWEEPS / "amplification_form" / "results" / "cells.csv",
            S.SWEEPS / "amplification_form" / "data" / "reference.csv",
            S.SWEEPS / "amplification_form" / "results" / "provenance.json",
            S.STRENGTHEN / "out" / "second" / "prediction.csv",
        ),
        "first": (
            S.STRENGTHEN / "raw" / "cells.csv",
            S.STRENGTHEN / "raw" / "reference.csv",
            S.SWEEPS / "product_amplification" / "results" / "provenance.json",
            S.STRENGTHEN / "out" / "first" / "prediction.csv",
        ),
    }
    modes = (
        ("none (0)", mpf(0)),
        ("const 2^-52", mpf(2) ** -52),
        ("const 2^-53", mpf(2) ** -53),
        ("ulp(value)", eta_ulp),
    )
    for name, paths in trees.items():
        cells_path, reference_path, provenance_path, recorded = paths
        states = S.recorded_states(recorded)
        print("=== %s experiment" % name)
        for label, eta in modes:
            model = S.Model(module, cells_path, reference_path, provenance_path, eta)
            inputs = model.band_inputs()
            winner = S.winners(model)
            aggregate, per = S.score(model, inputs, winner, 1.0)
            changed = []
            for band in ("A1", "A2"):
                for metric in ("absolute", "relative"):
                    rows = per[(band, metric)][3]
                    for n, ps, pd, sigma, state, measured, choice, right in rows:
                        if states[(band, metric, n)] != (state, choice):
                            changed.append("%s/%s/%d:%s->%s"
                                           % (band, metric, n, states[(band, metric, n)][0], state))
            print("  %-12s commits %3d right %3d abstains %3d; differs at %d: %s"
                  % (label, aggregate["commits"], aggregate["right"],
                     aggregate["abstains"], len(changed), ", ".join(changed)))
        print()


if __name__ == "__main__":
    main()
