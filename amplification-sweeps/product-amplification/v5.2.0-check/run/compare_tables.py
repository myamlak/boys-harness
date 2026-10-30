#!/usr/bin/env python3
"""Compare the regenerated tables and figure data with the deposited ones.

Run from the paper directory after the scratch table chain has completed:

    python sweeps/v5.2.0-check/run/compare_tables.py

Reads only files under this paper tree: the deposited tables/ and figures/
beside the manuscript, and the regenerated copies under
sweeps/v5.2.0-check/run/.  Line endings are folded to LF before comparing.
"""

import sys
from pathlib import Path

PAPER = Path(__file__).resolve().parents[3]
RUN = Path(__file__).resolve().parent


def folded(path):
    return path.read_bytes().replace(b"\r\n", b"\n").replace(b"\r", b"\n")


def compare(deposited, regenerated, label):
    a = folded(deposited)
    b = folded(regenerated)
    verdict = "SAME" if a == b else "DIFF"
    print("%-6s %-40s %8d bytes deposited" % (verdict, label, len(a)))
    if a != b:
        for index, (x, y) in enumerate(zip(a.split(b"\n"), b.split(b"\n")), 1):
            if x != y:
                print("   line %d" % index)
                print("     deposited: %r" % x[:200])
                print("     regenerated: %r" % y[:200])
                break
        if len(a.split(b"\n")) != len(b.split(b"\n")):
            print("   line counts differ: %d vs %d"
                  % (len(a.split(b"\n")) - 1, len(b.split(b"\n")) - 1))
    return a == b


def main():
    print("Tables and figure data: regenerated in the scratch tree against the")
    print("deposited tables beside the manuscript. CRLF folded to LF.")
    print("")
    ok = True
    for deposited in sorted((PAPER / "tables").iterdir()):
        regenerated = RUN / "tables" / deposited.name
        if not regenerated.is_file():
            print("MISSING regenerated %s" % deposited.name)
            ok = False
            continue
        ok &= compare(deposited, regenerated, "tables/" + deposited.name)
    deposited = PAPER / "figures" / "per-order-ratio.dat"
    regenerated = RUN / "figures" / "per-order-ratio.dat"
    if regenerated.is_file():
        ok &= compare(deposited, regenerated, "figures/per-order-ratio.dat")
    else:
        print("MISSING regenerated figures/per-order-ratio.dat")
        ok = False
    print("")
    print("verdict: %s" % ("all identical" if ok else "differences reported above"))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
