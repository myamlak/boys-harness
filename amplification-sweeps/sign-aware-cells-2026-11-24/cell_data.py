"""Load the four entries of the two chosen reversal cells.

Cell 1: A1 absolute order 2   (both delivered errors negative)
Cell 2: A2 absolute order 23  (mixed sign: chosen positive, other negative)

For each cell the chosen entry is the direct/product lane at x_prod and the
other entry is the recursed/ship lane at x_ship (harness_cells.csv, state
column).  Everything deposited is read; nothing is written outside this
directory.

Run:  python -I cell_data.py
"""
import csv
import math
import os
import struct
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
RECORD = REPO / "amplification-sweeps" / "certified-bound"
SWEEP = Path(os.environ.get("CERT_SWEEP", REPO / "amplification-sweeps" / "amplification-form"))
KERNEL = Path(os.environ.get("BOYS_ROOT", REPO.parent / "boys"))

sys.path.insert(0, str(RECORD / "code"))
import model                      # noqa: E402
import exactroute as er           # noqa: E402


CELLS = [(("A1", "absolute", 2)), (("A2", "absolute", 23))]


def bits(v):
    return struct.unpack(">Q", struct.pack(">d", v))[0]


def hexf(v):
    return float.hex(v)


def load_csv(path):
    with open(path, newline="") as h:
        return list(csv.DictReader(h))


def main():
    model.configure(kernel=KERNEL, exp_bits=RECORD / "exp" / "exp_bits.txt")
    harness = {((r["band"], r["metric"], int(r["n"]))): r
               for r in load_csv(RECORD / "results" / "harness_cells.csv")}
    cert = {(r["band"], float(r["x"]), int(r["n"]), r["method"]): r
            for r in load_csv(RECORD / "results" / "cells_cert.csv")}
    cells = {}
    for r in load_csv(SWEEP / "results" / "cells.csv"):
        cells[(r["band"], float(r["x"]), int(r["n"]))] = (
            float(r["recursed"]), float(r["direct"]))
    exp = model.load_exp()
    th = model.table()["thresholds"]

    def served(x):
        s = 0
        while s < 32 and x >= th[s + 1]:
            s += 1
        return s

    print("cell data")
    print("=" * 78)
    for key in CELLS:
        band, metric, n = key
        h = harness[key]
        xs, xp = float(h["x_ship"]), float(h["x_prod"])
        bi = 0 if band == "A1" else 1
        es = exp[xs]
        ds = cells[(band, xs, n)][0]
        dp = cells[(band, xp, n)][1]
        ms = model.ship_batch(xs, es[1])[n]
        mp_ = model.product_batch(xp, bi)[n]
        print(f"\n{band} {metric} n={n}  state={h['state']}  "
              f"loser={h['loser']}  trunc_order={h['trunc_order']}")
        print(f"  x_ship = {hexf(xs)}  ({xs!r})  served={served(xs)}  "
              f"mode={ms['mode']}")
        print(f"  x_prod = {hexf(xp)}  ({xp!r})  band_index={bi}  "
              f"mode={mp_['mode']}")
        print(f"  delivered ship bits = {bits(ds):016x}  {hexf(ds)}")
        print(f"  delivered prod bits = {bits(dp):016x}  {hexf(dp)}")
        print(f"  model ship bits     = {bits(ms['value']):016x}  "
              f"match={bits(ms['value']) == bits(ds)}")
        print(f"  model prod bits     = {bits(mp_['value']):016x}  "
              f"match={bits(mp_['value']) == bits(dp)}")
        print(f"  exp probe exp(-x)   = {hexf(es[0])}  "
              f"half = {hexf(es[1])}  bits={bits(es[1]):016x}")
        print(f"  half*2 == exp       = {es[1] * 2.0 == es[0]}")
        print(f"  record t_ship={h['t_ship']}  t_prod={h['t_prod']}")
        print(f"  record C_ship={h['C_ship']}  C_prod={h['C_prod']}")
        print(f"  record r_ship={h['r_ship']}  r_prod={h['r_prod']}")
        cs = cert[(band, xs, n, "ship")]
        cp = cert[(band, xp, n, "product")]
        print(f"  cert row ship: t={cs['t']} C={cs['C']} e={cs['e']} "
              f"r={cs['r']} mode={cs['mode']} c_seed={cs['c_seed']} "
              f"c_rec={cs['c_rec']} c_exp={cs['c_exp']}")
        print(f"  cert row prod: t={cp['t']} C={cp['C']} e={cp['e']} "
              f"r={cp['r']} mode={cp['mode']}")
        # exact route at 70 digits (record definition)
        sv = er.ship_values(xs, math.exp(-xs) / 2.0)[n]
        pv = er.product_values(xp)[n]
        print(f"  exactroute ship (exact e^-x/2) = {mp.nstr(sv, 25)}")
        print(f"  exactroute prod               = {mp.nstr(pv, 25)}")


if __name__ == "__main__":
    from mpmath import mp
    mp.dps = 70
    main()
