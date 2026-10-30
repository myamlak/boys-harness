#!/usr/bin/env python3
"""Provenance for the held-out sweep, without invoking git.

The revision this sweep ran at is the one the second sweep's record was taken
from: boys v5.1.1, commit 834bdbb18a30ac8538b519436429f78d41f187ce.  That
record, ../amplification_form/results/provenance.json, was written at the exact
tag with a clean tree and carries a sha256 for every kernel and reference file
it read.  This script does not run git.  It proves the same thing by content:

  * it reads the second sweep's record and names its tag and commit;
  * it digests every file under the kernel's include/ directory (the build
    reads headers only), the three reference-chain files, and the harness's
    own sources, grid and records;
  * it compares each digest the second sweep recorded with the current file
    and fails if any differs, so "the same kernel revision" is a content
    identity rather than a claim;
  * it records the kernel's version.hpp constants (5.1.1) and the constants
    the analyser reads (kTierThresholds, the region-A edges) from the headers
    themselves, through the same parser the other sweeps use.

The output is one JSON file beside the results.  No git argument is executed
here or anywhere in this sweep's pipeline.

Usage:
    python provenance.py --boys-tree C:\\Users\\myaml\\dev\\_boys-v511 \\
        --flags "<compile line>" --out results/provenance.json
"""

import argparse
import datetime
import json
import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "common"))
import boysprovenance  # noqa: E402

DEFAULT_TREE = boysprovenance.DEFAULT_TREE
KERNEL_RECORD = Path(__file__).resolve().parents[1] / "amplification_form/results/provenance.json"


def kernel_files(tree):
    root = tree / "include/boys"
    return sorted(path for path in root.rglob("*") if path.is_file())


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--boys-tree", default=os.environ.get("BOYS_TREE", DEFAULT_TREE))
    parser.add_argument("--kernel-record", default=str(KERNEL_RECORD))
    parser.add_argument("--harness", default="sweep_heldout")
    parser.add_argument("--harness-dir", default=".")
    parser.add_argument("--flags", default="")
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    tree = Path(args.boys_tree).resolve()
    if not (tree / "include/boys/version.hpp").is_file():
        raise SystemExit("not a boys tree: %s" % tree)

    record_path = Path(args.kernel_record)
    if not record_path.is_file():
        raise SystemExit("the second sweep's provenance record is missing: %s" % record_path)
    record = json.loads(record_path.read_text(encoding="utf-8"))
    recorded = record["boys"]
    if not recorded.get("tag_at_head") or recorded.get("dirty"):
        raise SystemExit("the second sweep's record does not name a clean exact tag")

    digests = {}
    compared = []
    mismatched = []
    for path in kernel_files(tree):
        relative = path.relative_to(tree).as_posix()
        entry = boysprovenance.digest(path)
        digests["boys/" + relative] = entry
        previous = record["digests"].get("boys/" + relative)
        if previous is None:
            compared.append((relative, "not in the second sweep's record"))
        elif previous["sha256"] != entry["sha256"]:
            mismatched.append(relative)
        else:
            compared.append((relative, "match"))

    for relative in boysprovenance.REFERENCE_FILES:
        path = tree / relative
        if path.is_file():
            entry = boysprovenance.digest(path)
            digests["boys/" + relative] = entry
            previous = record["digests"].get("boys/" + relative)
            if previous is None:
                compared.append((relative, "not in the second sweep's record"))
            elif previous["sha256"] != entry["sha256"]:
                mismatched.append(relative)
            else:
                compared.append((relative, "match"))

    root = Path(args.harness_dir).resolve()
    for path in sorted(root.rglob("*")):
        if path.is_file() and path.suffix in {".py", ".cpp", ".hpp", ".ps1", ".md"}:
            digests["harness/" + path.relative_to(root).as_posix()] = boysprovenance.digest(path)

    if mismatched:
        raise SystemExit(
            "the kernel tree does not match the second sweep's record at: %s" % ", ".join(mismatched)
        )

    matched = sum(1 for _, state in compared if state == "match")
    version_header = boysprovenance.read(tree / "include/boys/version.hpp")
    version = version_header.split('kVersionString = "', 1)[1].split('"', 1)[0]
    payload = {
        "generated_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "harness": args.harness,
        "flags": args.flags,
        "boys": {
            "tree": str(tree),
            "tag": recorded["tag_at_head"],
            "commit": recorded["commit"],
            "revision_source": str(record_path.resolve()),
            "revision_source_sha256": boysprovenance.digest(record_path)["sha256"],
            "version_header": version,
            "git_used": False,
            "content_check": {
                "files_in_record": len(record["digests"]),
                "kernel_files_digested": len(kernel_files(tree)),
                "matched": matched,
                "not_in_record": [name for name, state in compared if state != "match"],
                "mismatched": mismatched,
            },
        },
        "constants": boysprovenance.parse_constants(tree),
        "build_defaults": boysprovenance.parse_build_defaults(tree),
        "pinned_policy": boysprovenance.parse_pinned_policy(),
        "compiler": boysprovenance.compiler_banner(),
        "machine": boysprovenance.machine(),
        "digests": digests,
    }

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2, sort_keys=True)
        handle.write("\n")
    print(
        "wrote %s: kernel v%s at %s, %d kernel files matched against %s"
        % (out, version, recorded["commit"][:12], matched, record_path.name)
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
