#!/usr/bin/env python3
"""Check the deposited recurrence_arithmetic against the validated model.

The extended function in sweeps/amplification_form/analyze.py is the one both
sweeps' scoring imports.  This script compares its total, at every sampled
model cell, with score_exp.Model's rounding term plus the exponential term at
eta_E = 2^-52; the latter reproduced all 132 recorded states at eta_E = 0 and
its closed form was checked against the kernel's product form
(score-exp.log).  Writes nothing.
"""

import sys
from mpmath import mp, mpf
import score_exp as S

mp.dps = 70

TREES = {
    "second": (
        S.SWEEPS / "amplification_form" / "results" / "cells.csv",
        S.SWEEPS / "amplification_form" / "data" / "reference.csv",
        S.SWEEPS / "amplification_form" / "results" / "provenance.json",
    ),
    "first": (
        S.STRENGTHEN / "raw" / "cells.csv",
        S.STRENGTHEN / "raw" / "reference.csv",
        S.SWEEPS / "product_amplification" / "results" / "provenance.json",
    ),
}

ORDERS = (1, 2, 3, 5, 9, 16, 24, 32)


def main():
    module = S.load_analysis()
    worst = mpf(0)
    where = None
    checked = 0
    for name, (cells_path, reference_path, provenance_path) in TREES.items():
        model = S.Model(module, cells_path, reference_path, provenance_path, mpf(2) ** -52)
        model.band_inputs()
        keys = sorted(model.tables)
        for index, (band, x) in enumerate(keys):
            if index % 37:
                continue
            refs = model.reference[(band, x)]
            for n in ORDERS:
                if n > 32:
                    continue
                rounding, expt = model.arithmetic(n, band, x)
                expected = rounding + expt
                got = module.recurrence_arithmetic(n, mpf(x), refs)
                scale = abs(expected) if expected != 0 else mpf(1)
                deviation = abs(got - expected) / scale
                checked += 1
                if deviation > worst:
                    worst, where = deviation, (name, band, n, repr(x))
    print("checked %d sampled cells; worst relative deviation %s at %s"
          % (checked, mp.nstr(worst, 5), where))
    if worst > mpf(10) ** -58:
        raise SystemExit("deposited function and validated model disagree")
    print("agreement is at the script's 60-digit working precision")


if __name__ == "__main__":
    main()
