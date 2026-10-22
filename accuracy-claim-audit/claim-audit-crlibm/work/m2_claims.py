"""M2: measure crlibm against the independent MPFR reference inside -- and
just outside -- the domains its manual claims.

Every query is a (function, input) pair; the reference decides all four
rounding-mode results at once, so each query checks four documented claims at
once (the *_rn/_ru/_rd/_rz entry points).
"""
import os
import sys
import math
import random
import struct
import json
from multiprocessing import Pool

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ref
import run_crlibm

DBL_MAX = float.fromhex("0x1.fffffffffffffp+1023")
MIN_SUB = float.fromhex("0x0.0000000000001p-1022")
MAX_SUB = float.fromhex("0x0.fffffffffffffp-1022")
MIN_NORM = float.fromhex("0x1.0p-1022")
PI = math.pi

rnd = random.Random(20260924)


def hexof(x):
    return struct.pack(">d", x).hex()


def rand_in(lo, hi, n, logspace=False):
    out = []
    for _ in range(n):
        if logspace:
            a, b = math.log(abs(lo) + 1e-320), math.log(abs(hi))
            v = math.exp(rnd.uniform(a, b))
        else:
            v = lo + rnd.random() * (hi - lo)
        out.append(v if rnd.random() < 0.5 else -v)
    return out


def from_bits(hi, lo):
    return struct.unpack(">d", struct.pack(">II", hi, lo))[0]


def neighbours(x, k=4):
    out = [x]
    for i in range(1, k + 1):
        out.append(math.nextafter(x, math.inf))
        out.append(math.nextafter(x, -math.inf))
        x = out[-1]
    return out


def gen(func, n):
    """Function-specific adversarial input sample."""
    xs = []
    edge = {
        "exp": [0.0, MIN_SUB, -MIN_SUB, 1.0, -1.0, -708.0, -709.0, -745.0,
                -745.2, -1024.0, 709.0, 710.0, 0.5, -0.5, 2 ** -30],
        "expm1": [0.0, MIN_SUB, -MIN_SUB, 2 ** -54, -2 ** -54, 2 ** -30,
                  -2 ** -30, -1.0, 1.0, 709.0, -709.0],
        "log": [1.0, MIN_SUB, MAX_SUB, MIN_NORM, DBL_MAX, 2.0, 0.5, 10.0,
                1e300, 1.0 + 2 ** -52, 1.0 - 2 ** -53, math.pi],
        "log2": [1.0, 2.0, 0.5, MIN_SUB, DBL_MAX, 1.0 + 2 ** -52],
        "log10": [1.0, 10.0, 1e100, 1e-100, MIN_SUB, DBL_MAX] +
                 [10.0 ** k for k in range(-300, 301, 7)],
        "log1p": [0.0, MIN_SUB, -MIN_SUB, 2 ** -54, -2 ** -54, 1.0, -0.5,
                  DBL_MAX, -1.0 + 2 ** -53],
        "sin": [0.0, MIN_SUB, -MIN_SUB, PI, -PI, PI / 2, -PI / 2, 1.0,
                0.5, 2 ** -17, 2.0 + 4675 / 8192, 4.0, 100.0, 1e10, 1e300],
        "cos": [0.0, PI, -PI, PI / 4, 2 ** -25, 2 ** -18, 2 ** -17,
                12867 / 8192, 3.0, 1e10, 1e300],
        "tan": [0.0, MIN_SUB, 2 ** -25, 2 ** -18, math.atan(2), PI / 4,
                PI / 2 - 1e-10, 1.0, 1e10],
        "asin": [0.0, MIN_SUB, -MIN_SUB, 1.0, -1.0, 0.5, -0.5, 1.0 - 2 ** -53,
                 2 ** -30],
        "acos": [1.0, -1.0, 0.0, 0.5, MIN_SUB, 1.0 - 2 ** -53],
        "atan": [0.0, MIN_SUB, 1.0, -1.0, DBL_MAX, 1e300, math.tan(PI / 4)],
        "sinh": [0.0, MIN_SUB, 2 ** -30, -2 ** -30, 1.0, 23.0, 40.0, 709.0,
                 710.0, 710.4758600739439, -710.4758600739439],
        "cosh": [0.0, MIN_SUB, 1.0, 23.0, 40.0, 709.0, 710.0,
                 710.4758600739439, -710.4758600739439],
        "sinpi": [0.0, MIN_SUB, -MIN_SUB, 0.5, -0.5, 1.0, -1.0, 2.0,
                  2 ** -57, 2 ** -55, 0.25, 1e10, 1e300],
        "cospi": [0.0, 0.5, -0.5, 1.0, MIN_SUB, 2 ** -24, 1.0 / 3.0, 1e10],
        "tanpi": [0.0, MIN_SUB, 2 ** -25, 2 ** -5, 0.25, 1.0 / 3.0,
                  0.5 - 2 ** -53, -2 ** -20, 1e10],
        "atanpi": [0.0, MIN_SUB, 1.0, -1.0, math.tan(2 ** -25 * PI),
                   math.tan(2 ** -5 * PI), DBL_MAX, 1e-8],
    }[func]
    for e in edge:
        xs.extend(neighbours(e, 2))
    spans = {
        "exp": [(-745.2, 0.0), (0.0, 709.8), (0.0, 2 ** -30), (-2 ** -30, 0)],
        "expm1": [(-1.0, 0.0), (0.0, 1.0), (0.0, 2 ** -30), (2 ** -54, 2 ** -30)],
        "log": [(MIN_SUB, 1.0), (1.0, DBL_MAX), (1.0, 1.0 + 2 ** -40)],
        "log2": [(MIN_SUB, 1.0), (1.0, DBL_MAX)],
        "log10": [(MIN_SUB, 1.0), (1.0, DBL_MAX)],
        "log1p": [(-1.0, 0.0), (0.0, 1.0), (-2 ** -54, 2 ** -54)],
        "sin": [(0.0, PI), (PI, 100.0), (100.0, 1e12)],
        "cos": [(0.0, PI), (PI, 100.0), (100.0, 1e12)],
        "tan": [(0.0, 1.5), (1.5, 100.0)],
        "asin": [(-1.0, 1.0), (-0.5, 0.5)],
        "acos": [(-1.0, 1.0)],
        "atan": [(-1e3, 1e3), (-1.0, 1.0)],
        "sinh": [(0.0, 1.0), (1.0, 710.5), (-1.0, 0.0)],
        "cosh": [(0.0, 1.0), (1.0, 710.5)],
        "sinpi": [(0.0, 1.0), (1.0, 2.0), (0.0, 2 ** -50)],
        "cospi": [(0.0, 1.0), (0.0, 0.5)],
        "tanpi": [(0.0, 0.5), (2 ** -25, 2 ** -5), (2 ** -30, 2 ** -20)],
        "atanpi": [(-1.0, 1.0), (0.0, 10 * math.tan(2 ** -5 * PI))],
    }[func]
    per = max(1, n // max(1, len(spans)))
    for (lo, hi) in spans:
        xs.extend(rand_in(lo, hi, per, logspace=False))
    return [x for x in xs if not (x != x)]


def _decide(args):
    func, xhex = args
    rn, ru, rd, rz, st = ref.reference(func, xhex)
    return (func, xhex, rn, ru, rd, rz, st)


def run(func, nqueries, pool):
    xs = gen(func, nqueries)
    queries = [(func, hexof(x)) for x in xs]
    got = run_crlibm.evaluate([(func, "rn", h) for _, h in queries])
    got_ru = run_crlibm.evaluate([(func, "ru", h) for _, h in queries])
    got_rd = run_crlibm.evaluate([(func, "rd", h) for _, h in queries])
    got_rz = run_crlibm.evaluate([(func, "rz", h) for _, h in queries])
    decided = undecided = 0
    fails = []
    if pool is None:
        refs = [_decide(q) for q in queries]
    else:
        refs = pool.map(_decide, queries, chunksize=64)
    for i, (f, h) in enumerate(queries):
        _, _, rn, ru, rd, rz, st = refs[i]
        if st != "decided":
            undecided += 1
            continue
        decided += 1
        for mode, want, g in (("rn", rn, got[i]), ("ru", ru, got_ru[i]),
                              ("rd", rd, got_rd[i]), ("rz", rz, got_rz[i])):
            if g != want:
                fails.append((f, mode, h, want, g))
    return dict(func=func, queries=len(queries), decided=decided,
                undecided=undecided, fails=fails)


def main():
    funcs = sys.argv[1:] or ["exp", "log", "log2", "log10", "log1p", "expm1",
                             "sin", "cos", "tan", "asin", "acos", "atan",
                             "sinh", "cosh", "sinpi", "cospi", "tanpi",
                             "atanpi"]
    n = int(os.environ.get("MQ", "6000"))
    allfails = []
    summary = []
    with Pool(int(os.environ.get("MP", "8"))) as pool:
        for func in funcs:
            r = run(func, n, pool)
            print("%-8s queries=%-6d decided=%-6d undecided=%-5d failures=%d"
                  % (r["func"], r["queries"], r["decided"], r["undecided"],
                     len(r["fails"])), flush=True)
            for f in r["fails"][:6]:
                print("    FAIL %s %s x=%s want=%s got=%s" % f, flush=True)
            summary.append((r["func"], r["queries"], r["decided"],
                            r["undecided"], len(r["fails"])))
            allfails.extend(r["fails"])
    with open("m2_results.json", "w") as fh:
        json.dump({"summary": summary, "fails": allfails}, fh, indent=0)
    print("total failures:", len(allfails))


if __name__ == "__main__":
    main()
