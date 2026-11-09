#!/usr/bin/env python3
"""Compare the extracted boys v5.3.1 tag tree with the evaluated-file list.

usage: python compare_identity.py <tag-tree> <digests.txt> [<v5.1.1-tree>]

Without the optional tree it reports the match against the three digest
columns of the pre-tag v531-check record. With it, it reports the direct
byte-for-byte comparison against the measured v5.1.1 tag tree, which is the
comparison the 45/10/55 identity counts state.
"""

import hashlib
import pathlib
import sys

tag = pathlib.Path(sys.argv[1])
dig = pathlib.Path(sys.argv[2])
v511 = pathlib.Path(sys.argv[3]) if len(sys.argv) > 3 else None

rows = []
for line in dig.read_text().splitlines():
    if not line or line.startswith("#"):
        continue
    parts = line.split("\t")
    if len(parts) == 4:
        rows.append(parts)


def h(p):
    return hashlib.sha256(pathlib.Path(p).read_bytes()).hexdigest()


groups = {
    "headers": [r for r in rows if r[0].startswith("include/boys/")],
    "src": [r for r in rows if r[0].startswith("src/")],
    "all": rows,
}

print("evaluated files: %d (62 include/boys headers, 14 src)" % len(rows))
for name, group in groups.items():
    same = {
        "measured column": sum(1 for p, *ds in group if h(tag / p) == ds[0]),
        "v5.3.0 column": sum(1 for p, *ds in group if h(tag / p) == ds[1]),
        "v5.3.1 column": sum(1 for p, *ds in group if h(tag / p) == ds[2]),
    }
    if v511 is not None:
        same["v5.1.1 tag tree"] = sum(
            1 for p, *_ in group if h(tag / p) == h(v511 / p))
    print("%-8s %d files: %s" % (
        name, len(group), ", ".join("%s %d" % kv for kv in same.items())))

if v511 is not None:
    diff = [p for p, *_ in rows if h(tag / p) != h(v511 / p)]
    print("differing from the v5.1.1 tag tree: %d files" % len(diff))
    print("  " + "\n  ".join(diff))
