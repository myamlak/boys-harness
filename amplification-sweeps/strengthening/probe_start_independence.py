"""Diagnose the worst measurable start-independence cell.

    python probe_start_independence.py > probe_start_independence.log 2>&1

The start-independence check compares (F_n - ref_n) / (F_start - ref_start)
with A_n(x)/A_start(x).  That is a seed-only error model: it assumes the
delivered error at order n is the start's error carried by the recurrence.
The recurrence also injects rounding at every step, and the paper charges it
by sum_l 3*2^-53*F_l*A_n/A_l.  This script reads the deposited chain for the
worst measurable cell of check (2), A2 start=2 degree=20 n=26
x=10.910756753077845, and puts the measured numerator beside the seed term
and the accumulated-arithmetic charge.  If the charge is of the numerator's
size, the check's deviation is the breakdown of a seed-only model, not a
failure of the kernel's error accounting.
"""

import csv
import math
from pathlib import Path

from mpmath import mp, mpf

HERE = Path(__file__).resolve().parent
DEPOSIT = HERE.parent
BAND = "A2"
X_TEXT = "10.910756753077845"
START = 2
DEGREE = 20
N = 26
NMAX = 32
EPS = 2.220446049250313e-16


def amplification(n, x):
    if n == 0:
        return mpf(1)
    odd = math.factorial(2 * n) // (2**n * math.factorial(n))
    return mpf(odd) / (2 * mpf(x)) ** n


def read_reference(path):
    reference = {}
    with open(path, newline="") as handle:
        reader = csv.reader(handle)
        next(reader)
        for row in reader:
            reference[(row[0], float(row[1]))] = [mpf(v) for v in row[2:]]
    return reference


def main():
    mp.dps = 70
    reference = read_reference(DEPOSIT / "amplification-form" / "data" / "reference.csv")
    refs = reference[(BAND, float(X_TEXT))]

    chain = {}
    seed_row = None
    with open(HERE / "inputs" / "paths.csv", newline="") as handle:
        reader = csv.reader(handle)
        header = next(reader)
        for row in reader:
            if row[1] != BAND or row[2] != X_TEXT:
                continue
            if row[3] != str(START) or row[4] != str(DEGREE):
                continue
            if row[0] == "seed":
                seed_row = row
            else:
                chain[int(row[5])] = float(row[6].replace("nan(ind)", "nan"))
    if seed_row is not None:
        seed = float(seed_row[6].replace("nan(ind)", "nan"))
    else:
        # Seed-route rows are emitted only at x == 0; the start row carries
        # the seed for every other argument.
        seed = chain[START]
    print("header:", header)
    print("cell:", BAND, "x =", X_TEXT, "start =", START, "degree =", DEGREE)
    print("values read from the chain:", len(chain), "orders", min(chain), "..", max(chain))
    if seed != chain[START]:
        print("seed row and start row differ; using the seed row")
    chain[START] = seed

    x = mpf(X_TEXT)
    seed_error = mpf(seed) - refs[START]
    numerator = mpf(chain[N]) - refs[N]
    expected = amplification(N, x) / amplification(START, x)
    transfer = numerator / seed_error
    deviation = abs(abs(transfer) - expected) / expected
    seed_term = abs(seed_error) * expected
    charge = mpf(0)
    pieces = []
    for l in range(START + 1, N + 1):
        term = 3 * mpf(2) ** -53 * abs(mpf(chain[l])) * amplification(N, x) / amplification(l, x)
        charge += term
        pieces.append((l, float(term)))
    print()
    print("seed F_%d = %.17g, reference = %s" % (START, seed, mp.nstr(refs[START], 25)))
    print("F_%d = %.17g, reference = %s" % (N, chain[N], mp.nstr(refs[N], 25)))
    print("seed error (F_start - ref)      = %s" % mp.nstr(seed_error, 8))
    print("numerator (F_n - ref_n)         = %s" % mp.nstr(numerator, 8))
    print("closed-form A_n/A_start         = %s" % mp.nstr(expected, 8))
    print("measured transfer               = %s" % mp.nstr(transfer, 8))
    print("deviation                       = %.6g" % float(deviation))
    print()
    print("seed term |seed_error| * A_n/A_start = %s" % mp.nstr(seed_term, 8))
    print("accumulated-arithmetic charge       = %s" % mp.nstr(charge, 8))
    print("delivered numerator / seed term    = %s" % mp.nstr(abs(numerator) / seed_term, 8))
    print("delivered numerator / charge       = %s" % mp.nstr(abs(numerator) / charge, 8))
    print("charge / seed term                 = %s" % mp.nstr(charge / seed_term, 8))
    print()
    biggest = sorted(pieces, key=lambda item: item[1], reverse=True)[:5]
    print("the five largest charge terms (l, 3*2^-53*|F_l|*A_n/A_l):")
    for l, term in biggest:
        print("  l=%2d  %.6g" % (l, term))
    print()
    print("interpretation: the deposited check's template is the seed term alone.")
    print("The delivered numerator is that term plus the injected rounding; the")
    print("charge is the size the paper's own account predicts for the injected")
    print("part, and the two differ by the factor printed above.")


if __name__ == "__main__":
    main()
