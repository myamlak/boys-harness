#!/usr/bin/env python3
"""The execution check's sweep trees against the paper tree's records.

Every file the run reads from a sweep tree (harness scripts, the sweep
translation unit, grids, fits, 60-digit references and their validation
logs) is compared with the paper tree's copy, line endings folded.  The
results directory is the run's own output and is compared with the
deposited records by compare_sweeps.py instead; build products and caches
are not records and are excluded.
"""

import hashlib
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
CHECK = HERE.parent
PAPER = CHECK.parents[1]

TREES = ("product_amplification", "amplification_form")
SKIP_DIRS = {"build", "results", "__pycache__"}
SKIP_SUFFIXES = (".obj", ".exe", ".pdb", ".ilk", ".log")


def folded(path):
    data = path.read_bytes().replace(b"\r\n", b"\n").replace(b"\r", b"\n")
    return hashlib.sha256(data).hexdigest()


def files_of(root):
    return {
        path.relative_to(root).as_posix(): path
        for path in root.rglob("*")
        if path.is_file()
        and not any(part in SKIP_DIRS for part in path.relative_to(root).parts)
        and not path.name.endswith(SKIP_SUFFIXES)
    }


def main():
    out = Path(sys.argv[1]) if len(sys.argv) > 1 else None
    lines = [
        "The execution check's sweep trees against the paper tree's records",
        "Every input and harness file the run reads, CRLF folded. Build and",
        "results directories are excluded: the results are compared by",
        "compare_sweeps.py and the build products are not records.",
        "",
    ]
    same = differing = missing = extra = 0
    for tree in TREES:
        paper = files_of(PAPER / "sweeps" / tree)
        run = files_of(HERE / "sweeps" / tree)
        lines.append("=== " + tree)
        for rel in sorted(paper):
            if rel not in run:
                lines.append("MISSING  " + rel)
                missing += 1
                continue
            if folded(paper[rel]) == folded(run[rel]):
                lines.append("SAME     " + rel)
                same += 1
            else:
                lines.append("DIFF     " + rel)
                differing += 1
        for rel in sorted(set(run) - set(paper)):
            lines.append("EXTRA    " + rel)
            extra += 1
        lines.append("")
    lines.append("compared: %d identical, %d differing, %d missing, %d extra"
                 % (same, differing, missing, extra))
    text = "\n".join(lines) + "\n"
    if out is None:
        sys.stdout.write(text)
    else:
        out.write_text(text, encoding="utf-8", newline="\n")
    return 0 if differing == 0 and missing == 0 and extra == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
