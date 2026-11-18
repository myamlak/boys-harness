"""Compare the second-platform per_order.csv with the deposited one.

Usage: compare_per_order.py DEPOSITED SECOND_PLATFORM
Uses (band, n, method) as the key; compares at_x_abs and at_x_rel and the
worst values.  Prints the number of differing rows and every differing key.
"""
import csv
import sys


def load(path):
    rows = {}
    with open(path, newline="") as h:
        for r in csv.DictReader(h):
            rows[(r["band"], int(r["n"]), r["method"])] = r
    return rows


def main():
    a, b = load(sys.argv[1]), load(sys.argv[2])
    print("rows", len(a), len(b))
    diff_x = diff_abs = diff_rel = 0
    for key in sorted(a):
        ra, rb = a[key], b.get(key)
        if rb is None:
            print("missing", key)
            continue
        if ra["at_x_abs"] != rb["at_x_abs"]:
            diff_x += 1
            print("at_x_abs differs", key, ra["at_x_abs"], "->", rb["at_x_abs"],
                  "worst", ra["worst_abs"], "->", rb["worst_abs"])
        if ra["at_x_rel"] != rb["at_x_rel"]:
            diff_x += 1
            print("at_x_rel differs", key, ra["at_x_rel"], "->", rb["at_x_rel"],
                  "worst", ra["worst_rel"], "->", rb["worst_rel"])
        if ra["worst_abs"] != rb["worst_abs"]:
            diff_abs += 1
        if ra["worst_rel"] != rb["worst_rel"]:
            diff_rel += 1
    print("at_x differences", diff_x, "worst_abs text differences", diff_abs,
          "worst_rel text differences", diff_rel)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
