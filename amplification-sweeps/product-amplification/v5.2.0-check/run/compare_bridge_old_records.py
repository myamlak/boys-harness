#!/usr/bin/env python3
"""Establish the old side of the v5.2.0 bridge comparison.

The new side is the v5.2.0 run in _bridge-v520/runs; the old side is the
records in _bridge-v520/old-records.  This script shows where those came from
(byte identity with the v5.1.1 run's directory), compares them file by file
with the v5.1.1 files in the harness working tree, and compares them with the
v1.0.11 tag's chemistry-bridge/outputs, read from the repository's object
store without git.  The tag's set is reported as what it is: its own record
states version 5.0.0, and the tag has no -v511 files, so it is the v5.0.0 run,
not the old side of this comparison.

Line endings are folded (CRLF to LF) before text comparisons; .xz streams are
decompressed first.  A difference is printed with its first differing line.
"""

import hashlib
import lzma
import re
import sys
from datetime import datetime
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import harness_refs as h

COPY = Path(r"C:\Users\myaml\dev\_bridge-v520\old-records")
V511 = Path(r"C:\Users\myaml\dev\_bridge-v511\runs")
NEW = Path(r"C:\Users\myaml\dev\_bridge-v520\runs")
WORK = Path(r"C:\Users\myaml\dev\papers-harness-repo\chemistry-bridge\outputs")
TAG_DIR = "chemistry-bridge/outputs/"
VERSION = re.compile(r"^\s+version\s+(\S+)\s+\(boys::VersionString\)", re.M)


def fold(data):
    return data.replace(b"\r\n", b"\n")


def is_binary(name):
    return name.endswith(".bin")


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def version_of(name, data):
    match = VERSION.search(data.decode("utf-8", errors="replace"))
    return match.group(1) if match else "?"


def first_diff(name, old, new):
    if is_binary(name):
        return "binary: %d vs %d bytes" % (len(old), len(new))
    left = fold(old).decode("utf-8", errors="replace").splitlines()
    right = fold(new).decode("utf-8", errors="replace").splitlines()
    for number, (a, b) in enumerate(zip(left, right), 1):
        if a != b:
            return "line %d: %r vs %r" % (number, a[:120], b[:120])
    if len(left) != len(right):
        return "line count: %d vs %d" % (len(left), len(right))
    return "line endings only"


def main():
    out = Path(sys.argv[1]) if len(sys.argv) > 1 else None
    old = {p.name: p.read_bytes() for p in COPY.iterdir() if p.is_file()}
    names = sorted(old)

    lines = [
        "Old-side records for the v5.2.0 bridge comparison, "
        + datetime.now().astimezone().isoformat(timespec="seconds"),
        "",
        "1. The old side's provenance: %s against %s" % (COPY, V511),
    ]
    same = differ = missing = extra = 0
    extras = []
    for name in names:
        source = V511 / name
        if not source.is_file():
            missing += 1
            lines.append("  MISSING in v5.1.1 runs: %s" % name)
        elif source.read_bytes() == old[name]:
            same += 1
        else:
            differ += 1
            lines.append("  DIFFERS from v5.1.1 runs: %s" % name)
    for path in sorted(V511.iterdir()):
        if path.is_file() and path.name not in old:
            extra += 1
            extras.append(path.name)
            lines.append("  NOT COPIED from v5.1.1 runs: %s  %d bytes"
                         % (path.name, path.stat().st_size))
    lines.append("  %d files; %d byte-identical, %d differing, %d missing, "
                 "%d in the source and not copied."
                 % (len(names), same, differ, missing, extra))
    for name in extras:
        new_path = NEW / name
        if new_path.is_file():
            left, right = sha256(V511 / name), sha256(new_path)
            lines.append("  %s: v5.1.1 sha256 %s; v5.2.0 sha256 %s; %s"
                         % (name, left, right, "identical" if left == right
                            else "DIFFERENT"))
        else:
            lines.append("  %s: no v5.2.0 counterpart" % name)
    stated = sorted({version_of(n, old[n]) for n in names
                     if n.endswith(".boys_record.txt")})
    lines.append("  the old side's boys_record files state version %s"
                 % ", ".join(stated))
    lines.append("  (the two uncopied files are the cc-pVTZ dense tensors; "
                 "they are compared above by SHA-256, not by copy)")
    lines.append("")
    section1_ok = differ == 0 and missing == 0

    chain = h.resolve("refs/tags/v1.0.11")
    tag = h.walk(h.field(chain[-1][1], "tree"))
    tagged = {p[len(TAG_DIR):]: tag[p] for p in tag if p.startswith(TAG_DIR)}

    work = {p.name: p for p in WORK.iterdir() if p.is_file()} \
        if WORK.is_dir() else {}

    lines.append("2. Against the harness working tree, %s" % WORK)
    lines.append("   (the -v511 files are not in the v1.0.11 tree; they are "
                 "identified by path and hash, not by a tag)")
    counts = {"identical": 0, "via-plain": 0, "differs": 0, "missing": 0}
    for name in names:
        v511 = name.replace(".", "-v511.", 1)
        data = None
        via = v511
        if v511 in work:
            data = work[v511].read_bytes()
        elif (v511 + ".xz") in work:
            data = lzma.decompress(work[v511 + ".xz"].read_bytes())
            via = v511 + ".xz"
        elif name in work:
            data = work[name].read_bytes()
            via = name
        if data is None:
            counts["missing"] += 1
            lines.append("  MISSING %-46s (no %s)" % (name, v511))
            continue
        if fold(data) == fold(old[name]):
            if via == name:
                counts["via-plain"] += 1
                lines.append("  SAME via plain %-44s (no %s)" % (via, v511))
            else:
                counts["identical"] += 1
                lines.append("  SAME %-46s (%s)" % (name, via))
        else:
            counts["differs"] += 1
            lines.append("  DIFFERS %-43s %s" % (name, first_diff(name, old[name], data)))
    lines.append("  %d same by -v511 name, %d same via the plain file, "
                 "%d differing, %d missing."
                 % (counts["identical"], counts["via-plain"],
                    counts["differs"], counts["missing"]))
    work_stated = sorted({version_of(n, work[n].read_bytes())
                          for n in work if n.endswith("-v511.boys_record.txt")})
    lines.append("  the working tree's -v511 boys_record files state version %s"
                 % ", ".join(work_stated))
    lines.append("")

    lines.append("3. Against the v1.0.11 tag's %s" % TAG_DIR)
    lines.append("   (read from the object store at tag %s, commit %s)"
                 % (chain[0][0].split()[1],
                    chain[-1][0].split()[1] if len(chain) > 1 else "?"))
    tag_stated = sorted({version_of(n, h.read_object(tagged[n][1])[1])
                         for n in tagged if n.endswith(".boys_record.txt")})
    lines.append("   the tag's boys_record files state version %s"
                 % ", ".join(tag_stated))
    counts = {"identical": 0, "differs": 0, "missing": 0}
    cache = {}
    for name in names:
        entry = tagged.get(name) or tagged.get(name + ".xz")
        via = name if name in tagged else name + ".xz"
        if entry is None:
            counts["missing"] += 1
            lines.append("  MISSING %-46s (not deposited)" % name)
            continue
        if via not in cache:
            blob = h.read_object(entry[1])[1]
            cache[via] = lzma.decompress(blob) if via.endswith(".xz") else blob
        data = cache[via]
        if fold(data) == fold(old[name]):
            counts["identical"] += 1
            lines.append("  SAME %-46s (%s)" % (name, via))
        else:
            counts["differs"] += 1
            lines.append("  DIFFERS %-43s %s" % (name, first_diff(name, old[name], data)))
    lines.append("  %d same, %d differing, %d missing."
                 % (counts["identical"], counts["differs"], counts["missing"]))
    lines.append("")
    lines.append(
        "The tag's set is the v5.0.0 run: its own record states 5.0.0, and the "
        "tag carries no -v511 names.  The v5.1.1 run's records are not "
        "deposited in v1.0.11; only the v5.1.1-check comparison records are, "
        "under amplification-sweeps/product-amplification/v5.1.1-9353/.  The "
        "differences above are therefore between two library versions, "
        "reported as such, and the old side of the v5.2.0 comparison is the "
        "v5.1.1 run's own output, established by section 1.")

    text = "\n".join(lines) + "\n"
    if out is None:
        sys.stdout.write(text)
    else:
        out.write_text(text, encoding="utf-8", newline="\n")
    return 0 if section1_ok else 1


if __name__ == "__main__":
    sys.exit(main())
