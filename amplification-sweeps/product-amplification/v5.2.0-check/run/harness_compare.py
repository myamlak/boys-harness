#!/usr/bin/env python3
import hashlib
"""Compare the harness v1.0.10 and v1.0.11 trees, and v1.0.11 with the paper.

No git command is run: the trees are read from the repository's object store
by harness_refs.py.  The v1.0.11 bundle paths are mapped onto this paper
tree's sweep records and compared blob for blob.
"""

from pathlib import Path

import harness_refs as h

PAPER = Path(r"C:\Users\myaml\dev\papers-a-v5\papers\boys-crossover")

TAGS = {}
for tag in ("v1.0.10", "v1.0.11"):
    chain = h.resolve("refs/tags/" + tag)
    commit_label, commit_body = chain[-1]
    TAGS[tag] = (commit_label.split()[1], h.field(commit_body, "tree"))


def tree_files(tag):
    sha, _ = TAGS[tag]
    return h.walk(h.field(h.read_object(sha)[1], "tree"))


def blob_of(entry):
    return h.read_object(entry[1])[1]


def compare_trees():
    old = tree_files("v1.0.10")
    new = tree_files("v1.0.11")
    added = sorted(set(new) - set(old))
    removed = sorted(set(old) - set(new))
    changed = sorted(p for p in set(old) & set(new) if old[p][1] != new[p][1])
    print("v1.0.10 -> v1.0.11: %d -> %d files" % (len(old), len(new)))
    print("  added   %d: %s" % (len(added), added))
    print("  removed %d: %s" % (len(removed), removed))
    print("  changed %d:" % len(changed))
    for path in changed:
        a = blob_of(old[path])
        b = blob_of(new[path])
        print("    %-95s %8d -> %8d bytes" % (path, len(a), len(b)))
        if len(a) == len(b):
            diff = sum(1 for x, y in zip(a, b) if x != y)
            print("      same length; %d differing byte(s); CRLF counts %d -> %d"
                  % (diff, a.count(b"\r\n"), b.count(b"\r\n")))
    print()


def paper_path(bundle_path):
    prefix = "amplification-sweeps/product-amplification/v5.1.1-9353/"
    if not bundle_path.startswith(prefix):
        return None
    rest = bundle_path[len(prefix):]
    if rest.startswith("sweeps/"):
        return PAPER / rest
    if rest.startswith("v5.1.1-check/"):
        return PAPER / "sweeps" / rest
    return None


def compare_with_paper():
    files = tree_files("v1.0.11")
    prefix = "amplification-sweeps/product-amplification/v5.1.1-9353/"
    print("v1.0.11 bundle against this paper tree")
    same = diff = missing = 0
    unmapped = []
    for path in sorted(files):
        target = paper_path(path)
        if target is None:
            if path.startswith(prefix):
                unmapped.append(path[len(prefix):])
            continue
        blob = blob_of(files[path])
        if path.endswith(".xz"):
            target = Path(str(target)[:-3])
        if not target.is_file():
            print("  MISSING in paper: %s" % target)
            missing += 1
            continue
        if path.endswith(".xz"):
            import lzma
            blob = lzma.decompress(blob)
        paper_bytes = target.read_bytes()
        if blob == paper_bytes:
            same += 1
            continue
        diff += 1
        print("  DIFF %s" % path)
        print("    bundle %d bytes sha1 %s"
              % (len(blob), hashlib.sha1(b"blob %d\0" % len(blob) + blob).hexdigest()))
        print("    paper  %d bytes sha1 %s"
              % (len(paper_bytes),
                 hashlib.sha1(b"blob %d\0" % len(paper_bytes) + paper_bytes).hexdigest()))
        first = next((i for i, (x, y) in enumerate(zip(blob, paper_bytes)) if x != y), None)
        if first is not None:
            print("       first differing byte at %d: %r vs %r"
                  % (first, blob[first:first + 40], paper_bytes[first:first + 40]))
        if len(blob) != len(paper_bytes):
            print("       lengths differ; common prefix identical: %s"
                  % (blob[:min(len(blob), len(paper_bytes))]
                     == paper_bytes[:min(len(blob), len(paper_bytes))]))
    print("  compared: %d identical, %d differing, %d missing" % (same, diff, missing))
    print("  bundle entries with no paper-tree path (not compared): %d"
          % len(unmapped))
    for path in unmapped:
        print("    " + path)
    print()


def print_tag(tag):
    chain = h.resolve("refs/tags/" + tag)
    print("=== %s" % tag)
    for label, body in chain:
        print("  " + label)
        if label.startswith("tag"):
            print("    tag      %s" % h.field(body, "tag"))
            print("    object   %s" % h.field(body, "object"))
            print("    type     %s" % h.field(body, "type"))
        if label.startswith("commit"):
            tree = h.field(body, "tree")
            files = h.walk(tree)
            print("    tree     %s" % tree)
            print("    files    %d" % len(files))
            if "VERSION" in files:
                blob = h.read_object(files["VERSION"][1])[1]
                print("    VERSION  %s" % blob.decode(errors="replace").strip())
    print()


if __name__ == "__main__":
    for tag in ("v1.0.9", "v1.0.10", "v1.0.11"):
        print_tag(tag)
    compare_trees()
    compare_with_paper()
