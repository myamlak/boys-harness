"""Does the paper's own arithmetic charge explain the start-independence check?

    python diagnose_start_independence.py > diagnose_start_independence.log 2>&1

The check compares (F_n - ref_n) / (F_start - ref_start) with A_n/A_start.
That template is a seed-only error model: it assumes the delivered error at
order n is the start's error carried by the recurrence.  The recurrence also
rounds at every step; the paper charges that by sum_l 3*2^-53*|F_l|*A_n/A_l.

Write e_s for the start error, a = A_n/A_start, R for the rounding injected by
the steps after the start, and C for the charge.  The upward step is affine and
its amplification multiplies error by (2l-1)/(2x) per step, so the delivered
numerator is e_s*a + R and |R| <= C.  Hence

    |measured transfer - a| / a  =  |R| / |e_s a|  <=  C / |e_s a|.

The script evaluates that bound on every sampled cell of the check, using the
same sampling, the same reference and the same gate as check_transfers.py, so
its totals must match that log.  For a non-canonical start the charge runs over
the steps the chain actually takes, l = start+1 .. n, with the same 3*2^-53
constant the paper uses.  All error arithmetic is at 70 digits.
"""

import csv
import math
from pathlib import Path

from mpmath import mp, mpf

mp.dps = 70

HERE = Path(__file__).resolve().parent
INPUTS = HERE / "inputs"
DEPOSIT = HERE.parent
EPS = 2.220446049250313e-16
THETA0 = 1.0855252345349333
CHARGE_CONSTANT = 3 * mpf(2) ** -53


def amplification(n, x):
    if n == 0:
        return mpf(1)
    odd = math.factorial(2 * n) // (2**n * math.factorial(n))
    return mpf(odd) / (2 * x) ** n


def read_reference(path):
    reference = {}
    with open(path, newline="") as handle:
        reader = csv.reader(handle)
        next(reader)
        for row in reader:
            reference[(row[0], float(row[1]))] = [mpf(v) for v in row[2:]]
    return reference


def main():
    reference = read_reference(DEPOSIT / "amplification-form" / "data" / "reference.csv")
    amp_cache = {}

    def amp(n, x):
        key = (n, x)
        if key not in amp_cache:
            amp_cache[key] = amplification(n, x)
        return amp_cache[key]

    rows_total = 0
    cells = []
    current_key = None
    argument_index = -1
    chains = {}

    with open(INPUTS / "paths.csv", newline="") as handle:
        reader = csv.reader(handle)
        next(reader)
        for row in reader:
            rows_total += 1
            route, band, x_text, start_s, degree_s, n_s, value_s, pert_s, seed_s, transfer_s = row
            x = float(x_text)
            start = int(start_s)
            degree = int(degree_s)
            n = int(n_s)
            key = (band, x)

            if key != current_key:
                current_key = key
                argument_index += 1
                chains = {}

            if route == "seed" or x == 0.0:
                continue

            value = float(value_s.replace("nan(ind)", "nan"))

            if start not in (0, 32) and argument_index % 16 == 0:
                if n == start:
                    chains[(start, degree)] = {
                        "seed_error": mpf(value) - reference[(band, x)][start],
                        "amp": mpf(1),
                        "charge": mpf(0),
                    }
                else:
                    chain = chains.get((start, degree))
                    if chain is not None:
                        refs = reference[(band, x)]
                        ratio = mpf(2 * n - 1) / (2 * mpf(x))
                        chain["amp"] *= ratio
                        chain["charge"] = (chain["charge"] * ratio
                                           + CHARGE_CONSTANT * abs(mpf(value)))
                        seed_error = chain["seed_error"]
                        if seed_error != 0 and mp.isfinite(chain["amp"]):
                            expected = amp(n, x) / amp(start, x)
                            if expected > 0:
                                numerator = mpf(value) - refs[n]
                                transfer = numerator / seed_error
                                if mp.isfinite(transfer):
                                    seed_term = abs(seed_error) * expected
                                    floor = 16 * mpf(EPS) * abs(mpf(value))
                                    measurable = seed_term > floor
                                    deviation = abs(abs(transfer) - expected) / expected
                                    bound = chain["charge"] / seed_term
                                    cells.append((float(deviation), float(bound),
                                                  bool(measurable), band, x_text,
                                                  start, degree, n,
                                                  float(transfer), float(expected),
                                                  x >= THETA0,
                                                  float(abs(numerator) / abs(mpf(value)))))

    measurable = [c for c in cells if c[2]]
    inside = [c for c in measurable if c[10]]
    below = [c for c in measurable if not c[10]]
    print("paths rows scanned:", rows_total)
    print("sampled cells:", len(cells), "  (check_transfers.log: 132496)")
    print("measurable cells:", len(measurable), "  (98 306)")
    print()
    print("split by the model's population (x >= theta0 = %.16f):" % THETA0)
    for label, subset in (("inside", inside), ("below", below)):
        bad = [c for c in subset if c[0] > c[1]]
        worst = max(subset, key=lambda c: c[0]) if subset else None
        print("  %-6s population: measurable %6d, bound violations %6d, worst deviation %.6g%s"
              % (label, len(subset), len(bad), worst[0] if worst else float("nan"),
                 (" at %s x=%s start=%d deg=%d n=%d"
                  % (worst[3], worst[4], worst[5], worst[6], worst[7])) if worst else ""))
    print()
    print("bound |measured - a|/a <= charge/seed_term, on every measurable cell:")
    for label, subset in (("all", measurable), ("inside", inside), ("below", below)):
        bad = [c for c in subset if c[0] > c[1]]
        ratios = sorted(c[0] / c[1] for c in subset if c[1] > 0)
        print("  %-6s: violations %5d of %5d; deviation/bound median %.3g  p90 %.3g"
              "  p99 %.3g  max %.6g"
              % (label, len(bad), len(subset), ratios[len(ratios) // 2],
                 ratios[int(0.9 * (len(ratios) - 1))],
                 ratios[int(0.99 * (len(ratios) - 1))], ratios[-1]))
    print()
    print("the two regimes the seed-only template recognises, inside the population:")
    for label, subset in (("inside", inside), ("below", below)):
        ok = [c for c in subset if c[1] <= 1.0]
        brk = [c for c in subset if c[1] > 1.0]
        print("  %-6s: seed term above the charge (bound <= 1): %5d cells, deviation > 1: %4d"
              % (label, len(ok), sum(1 for c in ok if c[0] > 1.0)))
        print("          charge above the seed term (bound  > 1): %5d cells, deviation > 1: %4d"
              % (len(brk), sum(1 for c in brk if c[0] > 1.0)))
    print()
    worst = max(measurable, key=lambda c: c[0])
    print("worst measurable deviation overall: %.6g at %s x=%s start=%d degree=%d n=%d"
          % (worst[0], worst[3], worst[4], worst[5], worst[6], worst[7]))
    print("  inside the population: %s" % ("yes" if worst[10] else "no"))
    print("  its bound, charge/seed_term: %.6g" % worst[1])
    print("  its measured transfer and A_n/A_start: %.6g  %.6g"
          % (worst[8], worst[9]))
    inside_worst = max(inside, key=lambda c: c[0])
    print("worst measurable deviation inside the population: %.6g at %s x=%s start=%d degree=%d n=%d"
          % (inside_worst[0], inside_worst[3], inside_worst[4], inside_worst[5],
             inside_worst[6], inside_worst[7]))
    print("  its bound, charge/seed_term: %.6g" % inside_worst[1])
    print("  measured transfer and A_n/A_start: %.6g  %.6g"
          % (inside_worst[8], inside_worst[9]))
    print("  delivered relative error |F_n - ref| / |F_n|: %.3g" % inside_worst[11])
    print("  (check_transfers.log reports 264.029 at A2 start=2 deg=20 n=26 x=10.910756753077845)")
    print()
    print("the five population cells that most exceed the charge bound, by deviation/bound:")
    for c in sorted(inside, key=lambda c: c[0] / c[1] if c[1] else float("inf"),
                    reverse=True)[:5]:
        print("  deviation %.6g  bound %.6g  ratio %.6g  delivered rel err %.3g"
              "  %s x=%s start=%d deg=%d n=%d"
              % (c[0], c[1], c[0] / c[1] if c[1] else float("inf"), c[11],
                 c[3], c[4], c[5], c[6], c[7]))


if __name__ == "__main__":
    main()
