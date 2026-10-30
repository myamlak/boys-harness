#!/usr/bin/env python3
"""Unified diff between the paper tree's provenance recorder and the copy
the execution check used.

The check ran without invoking git, so its copy of common/boysprovenance.py
takes the worktree identity from the environment, where run/verify_pub_tree.py
put the values it read from the repository's object store.  The rest of the
recorder, and common/reference.py, are unmodified; both diffs are printed.
"""

import difflib
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
PAPER = HERE.parents[2]

FILES = ("boysprovenance.py", "reference.py")


def main():
    out = Path(sys.argv[1]) if len(sys.argv) > 1 else None
    lines = [
        "The execution check's common scripts against the paper tree's copies",
        "A diff here is a difference in what ran; it is reported, not smoothed.",
        "",
    ]
    for name in FILES:
        left = PAPER / "sweeps/common" / name
        right = HERE / "sweeps/common" / name
        diff = list(difflib.unified_diff(
            left.read_text(encoding="utf-8").splitlines(),
            right.read_text(encoding="utf-8").splitlines(),
            fromfile="paper sweeps/common/" + name,
            tofile="check run/sweeps/common/" + name,
            lineterm=""))
        lines.append("=== " + name + ": " +
                     ("identical" if not diff else
                      "%d diff line(s)" % len(diff)))
        lines.extend(diff)
        lines.append("")
    text = "\n".join(lines) + "\n"
    if out is None:
        sys.stdout.write(text)
    else:
        out.write_text(text, encoding="utf-8", newline="\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
