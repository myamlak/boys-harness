#!/usr/bin/env python3
"""Print the exact overflow-threshold comparison behind the crlibm verdict.

Writes raw/crlibm_threshold_context.txt.

The named input is read as its exact binary64 value.  DBL_MAX and the
round-to-nearest overflow midpoint T_RN = 2^1024 - 2^970 are written as
exact mpf values, and the exact sinh is compared with both.  This is the
arithmetic behind the expected bits in raw/crlibm_sinh_reference.csv.
"""
from pathlib import Path
import struct

import mpmath as mp

ROOT = Path(__file__).resolve().parent.parent
RAW = ROOT / "raw"
DPS = 80
BASE = 0x408633CE8FB9F87E


def main() -> None:
    mp.mp.dps = DPS
    x = struct.unpack(">d", struct.pack(">Q", BASE))[0]
    xm = mp.mpf(x)
    v = mp.sinh(xm)
    dbl_max = mp.mpf(2) ** 1024 - mp.mpf(2) ** 971
    t_rn = mp.mpf(2) ** 1024 - mp.mpf(2) ** 970

    lines = [
        "crlibm sinh overflow-threshold context",
        f"mpmath {mp.__version__}, mp.dps = {DPS}",
        f"input bits = 0x{BASE:016x}, exact double x = {mp.nstr(xm, 40)}",
        f"DBL_MAX = 2^1024 - 2^971 = {mp.nstr(dbl_max, 40)}",
        f"T_RN    = 2^1024 - 2^970 = {mp.nstr(t_rn, 40)}",
        f"exact sinh(x)          = {mp.nstr(v, 45)}",
        f"sinh(x) > DBL_MAX      = {v > dbl_max}",
        f"sinh(x) >= T_RN        = {v >= t_rn}",
        f"sinh(x) / DBL_MAX - 1  = {mp.nstr(v / dbl_max - 1, 12)}",
        f"sinh(x) / T_RN - 1     = {mp.nstr(v / t_rn - 1, 12)}",
        f"asinh(DBL_MAX)         = {mp.nstr(mp.asinh(dbl_max), 40)}",
        f"asinh(T_RN)            = {mp.nstr(mp.asinh(t_rn), 40)}",
        f"asinh(exact sinh) - x  = {mp.nstr(mp.asinh(v) - xm, 6)}",
    ]
    text = "\n".join(lines) + "\n"
    (RAW / "crlibm_threshold_context.txt").write_text(text)
    print(text, end="")


if __name__ == "__main__":
    main()
