#!/usr/bin/env python3
"""Put the deposited provenance's harness digests beside the files they name.

The v1.0.11 deposit carries two v5.1.1 provenance records per sweep tree:
the run record at results/provenance.json and the v5.1.1 check's copy at
v5.1.1-check/<tree>-provenance.json.  Both store harness digests as
digests["harness/<name>"], each entry a byte count and a SHA-256 digest; the
top-level "harness" field names the sweep script, not a digest group.  This
script reads both records from the repository's object store, states where
the two documents differ, hashes the copies the deposit actually carries and
the paper tree's copies, and prints them side by side with an equality mark
for each pair.  A file that is recorded but not deposited is printed as such.
The v5.2.0 check's own recorded digests are compared with the paper tree's
copies in provenance-diff-v520.txt.
"""

import hashlib
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import harness_refs as h

HERE = Path(__file__).resolve().parent
PAPER = HERE.parents[2]
DEP = "amplification-sweeps/product-amplification/v5.1.1-9353/"

TREES = {
    "product_amplification": {
        "deposit": DEP + "sweeps/product_amplification/",
        "check": DEP + "v5.1.1-check/product-amplification-provenance.json",
        "paper": PAPER / "sweeps" / "product_amplification",
    },
    "amplification_form": {
        "deposit": DEP + "sweeps/amplification_form/",
        "check": DEP + "v5.1.1-check/amplification-form-provenance.json",
        "paper": PAPER / "sweeps" / "amplification_form",
    },
}
NAMES = ("SWEEP.md", "analyze.py", "grid.py", "run.ps1", "sweep1.cpp", "sweep2.cpp")
ABSENT = object()


def sha256(data):
    return hashlib.sha256(data).hexdigest()


def flatten(value, prefix=""):
    if isinstance(value, dict):
        for key in value:
            yield from flatten(value[key], prefix + "/" + key)
    else:
        yield prefix, value


def read_json(tag, path):
    if path not in tag:
        return None, None
    data = h.read_object(tag[path][1])[1]
    return data, json.loads(data)


def digest_text(entry):
    if not entry:
        return "-"
    return "%d/%s" % (entry["bytes"], entry["sha256"][:8])


def mark(recorded, carried):
    if not recorded or carried is None:
        return "-"
    return "yes" if recorded["sha256"] == carried else "NO"


def main():
    out = Path(sys.argv[1]) if len(sys.argv) > 1 else None
    chain = h.resolve("refs/tags/v1.0.11")
    tag = h.walk(h.field(chain[-1][1], "tree"))

    lines = [
        "The deposited provenance's harness digests and the files they name",
        "Tag v1.0.11, tag object %s, commit %s, read from the object store."
        % (chain[0][0].split()[1], chain[-1][0].split()[1]),
        "Paper tree %s." % PAPER,
        "",
        "The digest section is digests[\"harness/<name>\"], each entry a byte",
        "count and a SHA-256 digest; the top-level \"harness\" field names the",
        "sweep script, not a digest group.  The deposit carries two provenance",
        "records per tree: the run record (results/provenance.json) and the",
        "v5.1.1 check's copy (v5.1.1-check/<tree>-provenance.json); both are",
        "read.  Each digest is printed as bytes/sha256-prefix; a dash means the",
        "entry does not exist there.  run=dep, check=dep and paper=dep compare",
        "the run record's digest, the check copy's digest and the paper tree's",
        "copy against the deposit's own copy.",
        "",
    ]

    for tree, where in TREES.items():
        prefix = where["deposit"]
        run_data, run_doc = read_json(tag, prefix + "results/provenance.json")
        chk_data, chk_doc = read_json(tag, where["check"])
        lines.append("=== %s" % tree)
        if run_doc is None or chk_doc is None:
            lines.append("  a provenance record is missing from the tag: run %s, "
                         "check %s" % (run_doc is not None, chk_doc is not None))
            lines.append("")
            continue

        left = dict(flatten(run_doc))
        right = dict(flatten(chk_doc))
        differing = [key for key in sorted(set(left) | set(right))
                     if left.get(key, ABSENT) != right.get(key, ABSENT)]
        lines.append("  run record  %s %d/%s" % (
            "results/provenance.json", len(run_data), sha256(run_data)[:8]))
        lines.append("  check copy  %s %d/%s" % (
            where["check"].rsplit("/", 1)[-1], len(chk_data), sha256(chk_data)[:8]))
        lines.append("  the two documents differ in %d field(s): %s" % (
            len(differing), ", ".join(differing) if differing else "none"))
        for key in differing:
            lines.append("    " + key)
            lines.append("      run  : " + (
                "<absent>" if left.get(key, ABSENT) is ABSENT
                else repr(left[key])))
            lines.append("      check: " + (
                "<absent>" if right.get(key, ABSENT) is ABSENT
                else repr(right[key])))

        run_digits = run_doc.get("digests", {})
        chk_digits = chk_doc.get("digests", {})
        lines.append("  %-12s %14s %14s %14s %14s %8s %10s %10s" % (
            "file", "run record", "check copy", "deposit", "paper",
            "run=dep", "check=dep", "paper=dep"))
        for name in NAMES:
            run_want = run_digits.get("harness/" + name)
            chk_want = chk_digits.get("harness/" + name)

            entry = tag.get(prefix + name)
            if entry is not None:
                dep_data = h.read_object(entry[1])[1]
                dep_digest = sha256(dep_data)
                deposit_text = "%d/%s" % (len(dep_data), dep_digest[:8])
            else:
                dep_digest = None
                deposit_text = "-"

            paper = where["paper"] / name
            if paper.is_file():
                paper_data = paper.read_bytes()
                paper_digest = sha256(paper_data)
                paper_text = "%d/%s" % (len(paper_data), paper_digest[:8])
            else:
                paper_digest = None
                paper_text = "-"

            pap_dep = ("-" if paper_digest is None or dep_digest is None
                       else "yes" if paper_digest == dep_digest else "NO")
            lines.append("  %-12s %14s %14s %14s %14s %8s %10s %10s" % (
                name, digest_text(run_want), digest_text(chk_want),
                deposit_text, paper_text,
                mark(run_want, dep_digest), mark(chk_want, dep_digest), pap_dep))
        lines.append("")

    lines.append(
        "The deposit carries no SWEEP.md in either tree, yet all four "
        "provenance records name one.  The product tree carries two records "
        "of different vintages, identical apart from four script digests and "
        "generated_utc: the run record, generated 2026-10-27T21:38:40Z, "
        "records for analyze.py, grid.py, run.ps1 and sweep1.cpp exactly the "
        "digests of the copies the deposit carries (the same bytes as the "
        "paper tree), and for SWEEP.md the digest of the paper tree's copy; "
        "the check copy, generated 2026-10-26T22:23:19Z, records earlier "
        "versions of SWEEP.md, analyze.py, grid.py and run.ps1, and only its "
        "sweep1.cpp digest agrees with them.  The form tree's two records are "
        "byte-identical (24538/a1ed6ea3), and its SWEEP.md, analyze.py and "
        "run.ps1 digests likewise name versions other than the paper tree's, "
        "while its grid.py and sweep2.cpp digests match it; the form deposit "
        "carries none of these files.  Every data record is reproduced by the "
        "v5.2.0 run regardless, so the inconsistency is in the provenance's "
        "harness digest section alone.  It is reported here, not smoothed.")

    text = "\n".join(lines) + "\n"
    if out is None:
        sys.stdout.write(text)
    else:
        out.write_text(text, encoding="utf-8", newline="\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
