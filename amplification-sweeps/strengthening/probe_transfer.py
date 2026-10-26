"""Probe the worst exact-perturbation cell with perturbations above one ULP.

    python probe_transfer.py > probe_transfer.log 2>&1

The summary's worst cell is A2 x=11.896188857019089, start=2, degree=2, n=25:
the deposited one-ULP perturbation gives a measured transfer of 2.8125 against
the closed form A_25/A_2 = 0.4280, a relative deviation of 5.57.  This script
replays the deposited chain bit for bit with Python's correctly rounded
float64 and math.fma, checks the replay against the deposited base values and
perturbation deltas, and then grows the seed perturbation from 1 to 2^20 ULP.
If the deviation falls with the perturbation size, the 5.57 is the rounding
divergence of two one-ULP-different chains, not a failure of the amplification
law.  The cell table is extracted from inputs/paths.csv into
inputs/probe_cell.csv on the first run and reused after that.
"""

import csv
import math
import struct
from pathlib import Path

from mpmath import mp, mpf

HERE = Path(__file__).resolve().parent
INPUTS = HERE / "inputs"
CELL = INPUTS / "probe_cell.csv"
TARGET_X = "11.896188857019089"
BAND = "A2"
START = 2
DEGREE = 2
NMAX = 32
NPROBE = 25


def amplification(n, x):
    if n == 0:
        return mpf(1)
    odd = math.factorial(2 * n) // (2**n * math.factorial(n))
    return mpf(odd) / (2 * mpf(x)) ** n


def extract_cell():
    if CELL.exists():
        return
    seed_row = None
    chain = {}
    with open(INPUTS / "paths.csv", newline="") as handle:
        reader = csv.reader(handle)
        next(reader)
        for row in reader:
            if row[1] != BAND or row[2] != TARGET_X:
                continue
            if row[3] != str(START) or row[4] != str(DEGREE):
                continue
            if row[0] == "seed":
                seed_row = row
            else:
                chain[int(row[5])] = row
    if seed_row is None:
        # The sweep emits a "seed" route row only at x == 0; for every other
        # argument the seed is the n == start row of the chain, whose
        # pert_delta is seed_delta and whose transfer_ulp is 1.
        seed_row = chain[START]
    with open(CELL, "w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["route", "band", "x", "start", "degree", "n",
                         "value", "pert_delta", "seed_delta", "transfer_ulp"])
        writer.writerow(seed_row)
        for n in sorted(chain):
            writer.writerow(chain[n])


def read_cell():
    rows = list(csv.reader(open(CELL, newline="")))
    header, body = rows[0], rows[1:]
    seed = None
    chain = {}
    for row in body:
        record = dict(zip(header, row))
        if record["route"] == "seed":
            seed = float(record["value"])
        else:
            chain[int(record["n"])] = {
                "value": float(record["value"]),
                "pert_delta": float(record["pert_delta"]),
                "transfer": float(record["transfer_ulp"]),
            }
    if seed is None:
        seed = chain[START]["value"]
    return seed, chain


def refined_divide(a, divisor, reciprocal):
    quotient = a * reciprocal
    return math.fma(math.fma(-quotient, divisor, a), reciprocal, quotient)


def upward_step(l, f, x, reciprocal, half_exp):
    return refined_divide((l - 0.5) * f - half_exp, x, reciprocal)


def ulp_count_shift(seed, k):
    """seed + k ULP, by mantissa addition, exact for these values."""
    bits = struct.unpack("<Q", struct.pack("<d", seed))[0]
    return struct.unpack("<d", struct.pack("<Q", bits + k))[0]


def main():
    mp.dps = 70
    extract_cell()
    seed, chain = read_cell()
    x = float(TARGET_X)
    half_exp = 0.5 * math.exp(-x)
    reciprocal = 1.0 / x
    ulp = struct.unpack("<d", struct.pack("<Q", struct.unpack("<Q", struct.pack("<d", seed))[0] + 1))[0] - seed
    closed = amplification(NPROBE, x) / amplification(START, x)
    closed_float = float(closed)

    print("cell:", BAND, "x =", TARGET_X, "start =", START, "degree =", DEGREE)
    print("seed = %.17g  ULP = %.17g  closed form A_25/A_2 = %.17g"
          % (seed, ulp, closed_float))
    print()

    # ---- replay at one ULP ------------------------------------------------
    base = seed
    pert = ulp_count_shift(seed, 1)
    replay = {START: (base, pert - base)}
    for n in range(START + 1, NMAX + 1):
        base = upward_step(n, base, x, reciprocal, half_exp)
        pert = upward_step(n, pert, x, reciprocal, half_exp)
        replay[n] = (base, pert - base)

    base_bad = pert_bad = 0
    for n in range(START + 1, NMAX + 1):
        if replay[n][0] != chain[n]["value"]:
            base_bad += 1
        if replay[n][1] != chain[n]["pert_delta"]:
            pert_bad += 1
    print("replay against the deposited one-ULP chain, n = %d..%d:" % (START + 1, NMAX))
    print("  base values differing from the record: %d" % base_bad)
    print("  perturbation deltas differing from the record: %d" % pert_bad)
    print("  deposited transfer at n=%d: %.17g; replay: %.17g"
          % (NPROBE, chain[NPROBE]["transfer"], replay[NPROBE][1] / (ulp)))
    print()

    # ---- grow the perturbation -------------------------------------------
    print("k-ULP seed perturbation, measured transfer_t = (F_t(seed + k ulp) - F_t(seed)) / (k ulp)")
    print("k        transfer at n=25   deviation from closed form   max deviation over n=%d..%d"
          % (START + 1, NMAX))
    for k in (1, 2, 4, 16, 64, 256, 1024, 4096, 16384, 65536, 2**20, 2**24):
        perturbed_seed = ulp_count_shift(seed, k)
        base = seed
        pert = perturbed_seed
        delta = {}
        for n in range(START + 1, NMAX + 1):
            base = upward_step(n, base, x, reciprocal, half_exp)
            pert = upward_step(n, pert, x, reciprocal, half_exp)
            delta[n] = (pert - base) / (k * ulp)
        worst = max(abs(delta[n] - float(amplification(n, x) / amplification(START, x)))
                    / float(amplification(n, x) / amplification(START, x))
                    for n in range(START + 1, NMAX + 1))
        deviation = abs(delta[NPROBE] - closed_float) / closed_float
        print("%-8d %.17g   %.6g                  %.6g" % (k, delta[NPROBE], deviation, worst))

    print()
    print("closed form A_n/A_2 at n=25: %.17g" % closed_float)
    print("the one-ULP measured transfer 2.8125 is 6.57x the closed form;")
    print("2.8125 is the deposited perturbation delta 6.0986e-19 over the seed ULP 2.1684e-19,")
    print("i.e. the two chains differ by 90 ULP of F_25, not by the 13.7 ULP the closed form asks for.")


if __name__ == "__main__":
    main()
