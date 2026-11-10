#!/usr/bin/env python3
"""Build the v5.3.1 kernel record the held-out rerun's provenance reads.

The held-out provenance proves its kernel by content: every file under
``include/boys`` and the three reference-chain files is hashed and compared
with a record's digests.  The v5.1.1 run read the second sweep's record, which
names 18 headers and the three reference files (the files that sweep read).
This rerun reads the v5.3.1 digest column of the certified-bound
``v531-check`` record, the column the exact-tag replay proved against the tag
tree at harness v1.0.35, so all 62 headers and 14 source files are compared.
The three reference-chain digests are carried from the second sweep's record;
the builder asserts the v5.3.1 tree hashes to them, which is the statement
that the reference chain did not move between the two kernels.

The builder refuses to write the record unless every digest checks out in the
tree it is pointed at.

Usage:
    python kernel_record_v531.py --tree C:\\Users\\myaml\\dev\\_a-exact-tag-replay \\
        --v531-digests <harness>/amplification-sweeps/certified-bound/v531-check/digests.txt \\
        --previous <sweeps>/amplification_form/results/provenance.json \\
        --out results/kernel_record_v531.json
"""

import argparse
import hashlib
import json
from pathlib import Path

TAG_OBJECT = "d1e63f71c6f2501aecc6340e31216a42c2a390b9"
COMMIT = "2f83e40a6504368b7b2f13682e29b4f24ab2a40c"
TREE = "344df3ecc1c13122f59a8a98a78afbff3c34337b"
REFERENCE_FILES = [
    "tools/gen_boys_accuracy_gate_reference.py",
    "tests/data/boys_accuracy_gate_reference.csv",
    "tests/data/boys_reference.csv",
]


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tree", required=True)
    parser.add_argument("--v531-digests", required=True)
    parser.add_argument("--previous", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    tree = Path(args.tree).resolve()
    digests_path = Path(args.v531_digests).resolve()
    previous_path = Path(args.previous).resolve()

    version = (tree / "include/boys/version.hpp").read_text(encoding="utf-8")
    if 'kVersionString = "5.3.1"' not in version:
        raise SystemExit("not the v5.3.1 tree: %s" % tree)

    digests = {}
    rows = 0
    for line in digests_path.read_text(encoding="utf-8").splitlines():
        if not line or line.startswith("#"):
            continue
        parts = line.split("\t")
        if len(parts) != 4:
            raise SystemExit("unexpected digest row: %r" % line)
        relative, wanted = parts[0], parts[3]
        path = tree / relative
        if not path.is_file():
            raise SystemExit("missing in the tree: %s" % relative)
        got = sha256(path)
        if got != wanted:
            raise SystemExit("digest mismatch at %s" % relative)
        digests["boys/" + relative] = {"sha256": got}
        rows += 1
    if rows != 76:
        raise SystemExit("expected 76 evaluated files, read %d" % rows)

    previous = json.loads(previous_path.read_text(encoding="utf-8"))
    for relative in REFERENCE_FILES:
        wanted = previous["digests"]["boys/" + relative]["sha256"]
        got = sha256(tree / relative)
        if got != wanted:
            raise SystemExit("reference-chain file moved: %s" % relative)
        digests["boys/" + relative] = {"sha256": got}

    payload = {
        "boys": {
            "tag_at_head": "v5.3.1",
            "commit": COMMIT,
            "tag_object": TAG_OBJECT,
            "tree": TREE,
            "dirty": False,
            "revision_source": str(digests_path),
            "revision_source_sha256": sha256(digests_path),
            "note": (
                "the 76 evaluated files are the v5.3.1 column of the "
                "certified-bound v531-check record, the column the exact-tag "
                "replay proved against the tag tree; the three reference-chain "
                "digests are the second sweep's record's, asserted here against "
                "this tree, so the reference chain is unchanged from v5.1.1"
            ),
        },
        "digests": digests,
    }
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n",
                   encoding="utf-8")
    print("wrote %s: %d evaluated files + %d reference files, tag v5.3.1"
          % (out, rows, len(REFERENCE_FILES)))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
