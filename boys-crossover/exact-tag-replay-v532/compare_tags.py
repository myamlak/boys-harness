#!/usr/bin/env python3
"""Compare two boys tag trees file by file over the evaluated set.

usage: python compare_tags.py <from-tree> <to-tree> <digests.txt>

The evaluated set is read from the pre-tag record's digests.txt, which lists
the 62 include/boys headers and the 14 src files of the measured tree.
"""

import hashlib
import pathlib
import sys


def h(path):
    return hashlib.sha256(pathlib.Path(path).read_bytes()).hexdigest()


source = pathlib.Path(sys.argv[1])
target = pathlib.Path(sys.argv[2])
digests = pathlib.Path(sys.argv[3])

paths = [line.split("\t")[0] for line in digests.read_text().splitlines()
         if line and not line.startswith("#")]
differing = [p for p in paths if h(source / p) != h(target / p)]
print("evaluated files: %d" % len(paths))
print("identical: %d" % (len(paths) - len(differing)))
print("differing: %d" % len(differing))
for path in differing:
    print("  " + path)
