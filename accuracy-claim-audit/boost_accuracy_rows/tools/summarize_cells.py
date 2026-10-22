"""Per-row cell statistics the verdict line does not carry.

The gate reports the worst error over every cell, which on a row whose function underflows
includes cells where the comparison is degenerate: the double value is zero and the reference is
not, so the relative difference is infinite while the absolute error is negligible.  Those cells
are counted and not judged, exactly as the method says, but the worst *judged* value is the one
that says how close a met row came to its bound.
"""

import re
import subprocess
import sys

CELL = re.compile(r"cell x=(\S+)\s+observed=(\S+)\s+floor=(\S+)\s+scale=(\S+)\s+(\S+)")
ROW = re.compile(r"^ (\S+)\s+(\S+)\s+\S+\s+\[")


def main() -> int:
    out = subprocess.run([sys.argv[1], "--rows=" + sys.argv[2], "--cells"],
                         capture_output=True, text=True).stdout  # exit 1 is the gate failing, which is a result
    rows, current = {}, None
    for line in out.splitlines():
        m = ROW.match(line)
        if m:
            current = m.group(1)
            rows[current] = []
            continue
        m = CELL.search(line)
        if m and current:
            rows[current].append((float(m.group(1)), float(m.group(2)), float(m.group(4)), m.group(5) == "judged"))
    print(f"{'id':<4} {'cells':>6} {'judged':>7} {'degenerate':>11} {'worst judged':>13} {'worst cell':>12}")
    for rid, cells in rows.items():
        degen = [c for c in cells if c[2] == 0.0]
        judged = [c for c in cells if c[3]]
        worst_judged = max((c[1] for c in judged), default=0.0)
        worst = max((c[1] for c in cells), default=0.0)
        print(f"{rid:<4} {len(cells):>6} {len(judged):>7} {len(degen):>11} "
              f"{worst_judged:>13.4g} {worst:>12.4g}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())