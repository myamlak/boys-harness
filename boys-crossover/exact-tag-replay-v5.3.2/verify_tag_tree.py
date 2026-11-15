#!/usr/bin/env python3
"""Verify an extracted boys tag tree against `git ls-tree -r`.

usage: python -I verify_tag_tree.py <tree> <ls-tree-output>

Every blob listed by git is read from the extracted tree, hashed as a git
blob and compared with the recorded object id. Exit 0 when all match.
"""
import hashlib
import pathlib
import sys

root = pathlib.Path(sys.argv[1])
listing = pathlib.Path(sys.argv[2])

count = 0
bad = 0
for line in listing.read_text(encoding="utf-8").splitlines():
    meta, path = line.split("\t")
    mode, kind, sha = meta.split()
    if kind != "blob":
        continue
    data = (root / path).read_bytes()
    count += 1
    if hashlib.sha1(b"blob %d\0" % len(data) + data).hexdigest() != sha:
        bad += 1
        print("MISMATCH %s" % path)
print("files %d mismatched %d" % (count, bad))
sys.exit(1 if bad else 0)
