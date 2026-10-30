#!/usr/bin/env python3
"""Field-by-field difference between the v5.1.1-check provenance copies and
the v5.2.0 execution-check records.

The deposited side is sweeps/v5.1.1-check/<tree>-provenance.json, the copy
the v5.1.1 run deposited; the new side is
run/sweeps/<tree>/results/provenance.json.  Line endings are not involved
(the comparison is on the parsed JSON).  Writes to the path given as the
first argument, or to stdout when no path is given.
"""

import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
CHECK = HERE.parent
PAPER = CHECK.parents[1]

TREES = {
    "product_amplification": (
        PAPER / "sweeps/v5.1.1-check/product-amplification-provenance.json",
        HERE / "sweeps/product_amplification/results/provenance.json"),
    "amplification_form": (
        PAPER / "sweeps/v5.1.1-check/amplification-form-provenance.json",
        HERE / "sweeps/amplification_form/results/provenance.json"),
}

ABSENT = object()


def flatten(value, prefix=""):
    if isinstance(value, dict):
        for key in value:
            yield from flatten(value[key], prefix + "/" + key)
    else:
        yield prefix, value


def main():
    out = Path(sys.argv[1]) if len(sys.argv) > 1 else None
    lines = [
        "Field-by-field difference between the v5.1.1-check provenance copies",
        "and the v5.2.0 execution-check records.  The deposited side is",
        "sweeps/v5.1.1-check/<tree>-provenance.json, the copy the v5.1.1 run",
        "deposited; the new side is run/sweeps/<tree>/results/provenance.json.",
        "Line endings are not involved (the comparison is on the parsed JSON).",
        "",
    ]
    for tree, (deposited_path, new_path) in TREES.items():
        deposited = dict(flatten(
            json.loads(deposited_path.read_text(encoding="utf-8"))))
        new = dict(flatten(
            json.loads(new_path.read_text(encoding="utf-8"))))
        differing = []
        for key in sorted(set(deposited) | set(new)):
            left = deposited.get(key, ABSENT)
            right = new.get(key, ABSENT)
            if left is ABSENT or right is ABSENT or left != right:
                differing.append((key, left, right))
        lines.append("=== %s: %d differing field(s)" % (tree, len(differing)))
        for key, left, right in differing:
            lines.append("  " + key)
            lines.append("    deposited: " + (
                "<absent>" if left is ABSENT else repr(left)))
            lines.append("    new      : " + (
                "<absent>" if right is ABSENT else repr(right)))
        lines.append("")
    text = "\n".join(lines) + "\n"
    if out is None:
        sys.stdout.write(text)
    else:
        out.write_text(text, encoding="utf-8", newline="\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
