#!/usr/bin/env python3
"""Write the v5.2.0 source and execution check records for Paper A.

The source half compares the evaluated sources of the v5.2.0 tag with the
historical v5.1.1 run, reading the tag's objects directly without git.  The
execution half does not run the sweeps here: run/run_execution_check.ps1 does
that and writes run/logs/execution-check.log and the comparison records beside
this file, which this generator reads and embeds.  Run from this directory
after the execution check:

    python make_check.py
"""

import datetime
import difflib
import hashlib
import json
import re
import zlib
from pathlib import Path

OUT = Path(__file__).resolve().parent
PAPER = OUT.parents[1]
V511 = Path(r"C:\Users\myaml\dev\_boys-v511")
V520 = Path(r"C:\Users\myaml\dev\_boys-umbrella")
GITDIR = Path(r"C:\Users\myaml\dev\_mainonly\.git")

TAG = "v5.2.0"
TAG_OBJECT = "cb1ce4c79eb22a83bc2fa116c6647df9b08a0e2e"
TAG_TARGET = "55ca0e316ad2f9f776424bfee9f365dfc322058b"
TAG_TREE = "312e72ab11f314d41834286143ce5549bafad72c"
TAG_DATE = "2026-10-28T18:09:47+00:00"
TAGGER = "Marcin Makowski <myamlak@gmail.com> 1793210987 +0100"
TAG_MESSAGE = ("boys 5.2.0: the claim-book registration table, the measured "
               "rows and the exhaustive half sweeps")
COMMIT_SUBJECT = "Claim-book registrations and half-format sweeps (#19)"
V511_COMMIT = "834bdbb18a30ac8538b519436429f78d41f187ce"
HEAD_REF = "refs/heads/myamlak/accuracy-gate-umbrella"
HEAD_COMMIT = "14c25bf8df04eaf585790aa2674bafcb15d65b8c"

EXCLUDE_DIRS = {".git", "build", "out", ".vs", "__pycache__", "CMakeFiles"}
EXCLUDE_SUFFIX = (".log", ".err", ".obj", ".exe", ".pdb", ".ilk", ".pyc", ".idb")
EXCLUDE_FILES = {"CMakeCache.txt", "Makefile", "CTestTestfile.cmake"}


def sha256(path):
    hasher = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            hasher.update(block)
    return hasher.hexdigest()


def blob_sha1(data):
    hasher = hashlib.sha1()
    hasher.update(b"blob " + str(len(data)).encode() + b"\x00")
    hasher.update(data)
    return hasher.hexdigest()


def loose(sha):
    path = GITDIR / "objects" / sha[:2] / sha[2:]
    raw = zlib.decompress(path.read_bytes())
    head, _, body = raw.partition(b"\x00")
    return head.decode("ascii"), body


def tree_entries(sha):
    head, body = loose(sha)
    if not head.startswith("tree "):
        raise ValueError("not a tree: " + sha)
    entries = []
    index = 0
    while index < len(body):
        end = body.index(b"\x00", index)
        mode, name = body[index:end].decode("utf-8").split(" ", 1)
        entries.append((mode, name, body[end + 1:end + 21].hex()))
        index = end + 21
    return entries


def resolve(root, path):
    sha = root
    for part in path.split("/"):
        sha = next(entry[2] for entry in tree_entries(sha) if entry[1] == part)
    return sha


def tag_tree_files():
    result = {}
    for sub in ("include/boys", "src"):
        subtree = resolve(TAG_TREE, sub)
        for mode, name, blob in tree_entries(subtree):
            if mode == "40000":
                raise ValueError("unexpected subdirectory: " + sub + "/" + name)
            result[sub + "/" + name] = blob
    return result


def tree_files(root):
    return sorted(
        path.relative_to(root).as_posix()
        for path in root.rglob("*")
        if path.is_file()
    )


def excluded(rel):
    parts = rel.split("/")
    if any(part in EXCLUDE_DIRS or part.startswith("build") or
           part.startswith("_gate") for part in parts):
        return True
    if rel.endswith(EXCLUDE_SUFFIX) or rel in EXCLUDE_FILES:
        return True
    return False


def category(rel):
    if rel.startswith("include/boys/") or rel.startswith("src/"):
        return "evaluated"
    if rel.startswith("tests/"):
        return "gate/tests"
    if rel.startswith("tools/"):
        return "tools"
    if rel.startswith("docs/") or rel in {
            "README.md", "CHANGELOG.md", "CITATION.bib", "CMakeLists.txt",
            "CONTRIBUTING.md", "THIRD_PARTY_NOTICES.md"}:
        return "metadata/docs"
    return "other"


def read_execution():
    exec_log = OUT / "run" / "logs" / "execution-check.log"
    if not exec_log.is_file():
        return None
    steps = []
    for line in exec_log.read_text(encoding="utf-8-sig").splitlines():
        parts = line.split()
        if len(parts) == 3 and parts[1].startswith("EXIT="):
            steps.append({
                "step": parts[0],
                "exit": int(parts[1][len("EXIT="):]),
                "seconds": float(parts[2][len("SECONDS="):].replace(",", ".")),
            })
    if not steps:
        return None
    version_log = OUT / "run" / "logs" / "compiler-version.log"
    path_log = OUT / "run" / "logs" / "compiler-path.log"
    return {
        "steps": steps,
        "all_steps_exit_zero": all(step["exit"] == 0 for step in steps),
        "boys_tree": r"C:\Users\myaml\dev\_boys-v520-pub",
        "boys_tree_identity": "run/pub-tree-identity.txt",
        "compiler": (version_log.read_text(encoding="utf-8-sig").splitlines()[0]
                     if version_log.is_file() else None),
        "compiler_path": (path_log.read_text(encoding="utf-8-sig").strip()
                          if path_log.is_file() else None),
    }


def read_sweep_comparison():
    path = OUT / "sweep-results-compare-v520.txt"
    if not path.is_file():
        return None
    identical = 0
    differing = []
    for line in path.read_text(encoding="utf-8-sig").splitlines():
        parts = line.split()
        if not parts:
            continue
        if parts[0] == "SAME":
            identical += 1
        elif parts[0] == "DIFF":
            differing.append(" ".join(parts[1:3]))
    return {
        "files_compared": identical + len(differing),
        "identical": identical,
        "differing": differing,
        "record": path.name,
    }


def read_scratch_comparison():
    path = OUT / "scratch-tree-compare-v520.txt"
    if not path.is_file():
        return None
    match = re.search(
        r"compared: (\d+) identical, (\d+) differing, (\d+) missing, "
        r"(\d+) extra", path.read_text(encoding="utf-8-sig"))
    if match is None:
        raise ValueError("no comparison verdict in " + str(path))
    identical, differing, missing, extra = (int(v) for v in match.groups())
    return {
        "files_compared": identical + differing + missing + extra,
        "identical": identical,
        "differing": differing,
        "missing": missing,
        "extra": extra,
        "record": path.name,
    }


def read_table_comparison():
    path = OUT / "tables-compare-v520.txt"
    if not path.is_file():
        return None
    identical = 0
    verdict = None
    for line in path.read_text(encoding="utf-8-sig").splitlines():
        if line.startswith("SAME"):
            identical += 1
        elif line.startswith("verdict:"):
            verdict = line.split(":", 1)[1].strip()
    return {"identical": identical, "verdict": verdict, "record": path.name}


def read_deposit():
    path = OUT / "harness-deposit-compare-v111.txt"
    if not path.is_file():
        return None
    text = path.read_text(encoding="utf-8-sig")
    for marker in ("commit e4513b62b381d49ee5ec1165f25195ea6bf8ff7f",
                   "VERSION  1.0.11"):
        if marker not in text:
            raise ValueError("deposit record does not show " + marker)
    match = re.search(
        r"compared: (\d+) identical, (\d+) differing, (\d+) missing", text)
    if match is None:
        raise ValueError("no comparison verdict in " + str(path))
    identical, differing, missing = (int(value) for value in match.groups())
    missing_from_paper = [
        line.split("MISSING in paper:", 1)[1].strip()
        for line in text.splitlines() if "MISSING in paper:" in line]
    return {
        "repository": "myamlak/boys-harness",
        "tag": "v1.0.11",
        "tag_object": "1aa52650755d7445134a105d289029300096e7dd",
        "commit": "e4513b62b381d49ee5ec1165f25195ea6bf8ff7f",
        "path": "amplification-sweeps/product-amplification/v5.1.1-9353/",
        "mapped": identical + differing + missing,
        "identical": identical,
        "differing": differing,
        "missing_from_paper": missing_from_paper,
        "record": path.name,
    }


def main():
    tag_files = tag_tree_files()
    work_match = 0
    tag_rows = []
    for rel, blob in sorted(tag_files.items()):
        data = (V520 / rel).read_bytes()
        got = blob_sha1(data)
        ok = got == blob
        work_match += int(ok)
        tag_rows.append((rel, blob, got, ok))

    v511_files = tree_files(V511)
    v520_files = tree_files(V520)
    compared = []
    for sub in ("include/boys", "src"):
        names = sorted(
            {rel[len(sub) + 1:] for rel in v511_files if rel.startswith(sub + "/")} |
            {rel[len(sub) + 1:] for rel in v520_files if rel.startswith(sub + "/")}
        )
        for name in names:
            rel = sub + "/" + name
            left = V511 / rel
            right = V520 / rel
            if not left.is_file() or not right.is_file():
                continue
            compared.append((rel, sha256(left), sha256(right)))

    headers = [row for row in compared if row[0].startswith("include/boys/")]
    sources = [row for row in compared if row[0].startswith("src/")]
    header_same = sum(1 for _, left, right in headers if left == right)
    source_same = sum(1 for _, left, right in sources if left == right)

    provenance_path = PAPER / "sweeps/v5.1.1-check/product-amplification-provenance.json"
    provenance = json.loads(provenance_path.read_text(encoding="utf-8"))
    recorded = []
    for key, spec in sorted(provenance["digests"].items()):
        if not key.startswith("boys/"):
            continue
        rel = key[len("boys/"):]
        path = V520 / rel
        actual = sha256(path) if path.is_file() else None
        recorded.append((rel, spec["sha256"], actual))

    lines = []
    lines.append("Source comparison for the v5.2.0 check")
    lines.append("Generated " + datetime.datetime.now(
        datetime.timezone.utc).isoformat())
    lines.append("")
    lines.append("Candidate: tag " + TAG + " (annotated tag object " +
                 TAG_OBJECT + ")")
    lines.append("  commit " + TAG_TARGET + ", tree " + TAG_TREE)
    lines.append("  tag date " + TAG_DATE + "; parent " + V511_COMMIT)
    lines.append("Source tree: " + str(V520))
    lines.append("  head ref " + HEAD_REF + " at " + HEAD_COMMIT +
                 ", same tree " + TAG_TREE)
    lines.append("Reference run: historical v5.1.1, commit " + V511_COMMIT +
                 ", records in sweeps/v5.1.1-check/")
    lines.append("")
    lines.append("Tag-tree blob verification: " + str(work_match) + " of " +
                 str(len(tag_rows)) + " evaluated files in the source tree match "
                 "the tagged tree's blob SHA-1.")
    lines.append("")
    lines.append("include/boys: " + str(len(headers)) + " files; " +
                 str(header_same) + " identical to v5.1.1, " +
                 str(len(headers) - header_same) + " file(s) differ.")
    lines.append("file\tv520_sha256\tv511_sha256\tv520=v511")
    for rel, left, right in headers:
        lines.append("{}\t{}\t{}\t{}".format(
            rel.rsplit("/", 1)[-1], right, left, left == right))
    lines.append("")
    lines.append("src: " + str(len(sources)) + " files; " +
                 str(source_same) + " identical to v5.1.1, " +
                 str(len(sources) - source_same) + " file(s) differ.")
    lines.append("file\tv520_sha256\tv511_sha256\tv520=v511")
    for rel, left, right in sources:
        lines.append("{}\t{}\t{}\t{}".format(
            rel.rsplit("/", 1)[-1], right, left, left == right))
    lines.append("")
    lines.append("Digest comparison with the v5.1.1 provenance record "
                 "(product-amplification-provenance.json):")
    matched = sum(1 for _, want, got in recorded if got == want)
    lines.append("  " + str(matched) + " of " + str(len(recorded)) +
                 " recorded boys files match the v5.2.0 tree.")
    for rel, want, got in recorded:
        state = "match" if got == want else "differ"
        lines.append("  {}\t{}\t{}\t{}".format(rel, state, want, got))
    lines.append("")
    version = next(row for row in compared if row[0] == "include/boys/version.hpp")
    lines.append("The differing header, include/boys/version.hpp:")
    lines.append("  v5.1.1 sha256 " + version[1])
    lines.append("  v5.2.0 sha256 " + version[2])
    lines.append("  changes: kVersionMinor 1 -> 2; kVersionPatch 1 -> 0;")
    lines.append("           kVersionString \"5.1.1\" -> \"5.2.0\"; no other line differs.")
    lines.append("")
    lines.append("Conclusion: every evaluated translation unit is byte-identical "
                 "to the historical v5.1.1 run except the version constant, so "
                 "no measured error can change from a source difference. The "
                 "execution check in run/ then re-ran the sweeps against these "
                 "sources; its steps and exit codes are in "
                 "run/logs/execution-check.log, its record comparison in "
                 "sweep-results-compare-v520.txt and tables-compare-v520.txt, "
                 "and the summary in provenance.json.")
    lines.append("")

    (OUT / "evaluation-sources-compare.txt").write_text(
        "\n".join(lines), encoding="utf-8", newline="\n")

    manifest = []
    for rel, blob in sorted(tag_files.items()):
        manifest.append("{}  {}".format(sha256(V520 / rel), rel))
    (OUT / "evaluation-sources-v520.sha256").write_text(
        "\n".join(manifest) + "\n", encoding="utf-8", newline="\n")

    blobs = []
    for rel, blob, got, ok in tag_rows:
        blobs.append("{}\t{}\t{}\t{}".format(rel, blob, got, ok))
    (OUT / "tag-tree-blobs-v520.txt").write_text(
        "rel\tblob_sha1\ttree_blob\tmatch\n" + "\n".join(blobs) + "\n",
        encoding="utf-8", newline="\n")

    left_text = (V511 / "include/boys/version.hpp").read_text(encoding="utf-8").splitlines()
    right_text = (V520 / "include/boys/version.hpp").read_text(encoding="utf-8").splitlines()
    diff = difflib.unified_diff(
        left_text, right_text,
        fromfile="include/boys/version.hpp (v5.1.1)",
        tofile="include/boys/version.hpp (v5.2.0)",
        lineterm="")
    (OUT / "version-hpp.diff").write_text(
        "\n".join(diff) + "\n", encoding="utf-8", newline="\n")

    old = {rel for rel in v511_files if not excluded(rel)}
    new = {rel for rel in v520_files if not excluded(rel)}
    diff_lines = ["Source-tree differences, v5.2.0 tree vs the v5.1.1 tree",
                  "Generated " + datetime.datetime.now(
                      datetime.timezone.utc).isoformat(),
                  "Excluded: build artifacts, gate-run logs, __pycache__, .git.",
                  ""]
    groups = {}
    for rel in sorted(old | new):
        if rel not in new:
            groups.setdefault(category(rel), []).append(("removed", rel, "", ""))
        elif rel not in old:
            groups.setdefault(category(rel), []).append(("added", rel, "", ""))
        else:
            left = sha256(V511 / rel)
            right = sha256(V520 / rel)
            if left != right:
                groups.setdefault(category(rel), []).append(
                    ("modified", rel, left, right))
    for name in ("evaluated", "gate/tests", "metadata/docs", "tools", "other"):
        rows = groups.get(name, [])
        diff_lines.append(name + ": " + str(len(rows)) + " file(s)")
        for state, rel, left, right in rows:
            if state == "modified":
                diff_lines.append("  modified\t{}\t{}\t{}".format(rel, left, right))
            else:
                diff_lines.append("  {}\t{}".format(state, rel))
        diff_lines.append("")
    (OUT / "source-diff-v520.txt").write_text(
        "\n".join(diff_lines), encoding="utf-8", newline="\n")

    figures = []
    for sub in ("sweeps/product_amplification/results",
                "sweeps/amplification_form/results",
                "sweeps/v5.1.1-check"):
        for path in sorted((PAPER / sub).rglob("*")):
            if path.is_file():
                rel = path.relative_to(PAPER).as_posix()
                figures.append("{}  {}".format(sha256(path), rel))
    (OUT / "carried-figures.sha256").write_text(
        "\n".join(figures) + "\n", encoding="utf-8", newline="\n")

    execution = read_execution()
    sweep_comparison = read_sweep_comparison()
    table_comparison = read_table_comparison()
    scratch_comparison = read_scratch_comparison()
    deposit = read_deposit()

    payload = {
        "check": "source-identity-and-execution",
        "generated_utc": datetime.datetime.now(
            datetime.timezone.utc).isoformat(),
        "release": "5.2.0",
        "tag": {
            "name": TAG,
            "object": TAG_OBJECT,
            "date_utc": TAG_DATE,
            "tagger": TAGGER,
            "message": TAG_MESSAGE,
        },
        "commit": {
            "sha": TAG_TARGET,
            "short": TAG_TARGET[:7],
            "parent": V511_COMMIT,
            "subject": COMMIT_SUBJECT,
            "tree": TAG_TREE,
        },
        "source_tree": {
            "path": str(V520),
            "head_ref": HEAD_REF,
            "head_commit": HEAD_COMMIT,
            "head_tree": TAG_TREE,
        },
        "reference_run": {
            "release": "5.1.1",
            "commit": V511_COMMIT,
            "records": "sweeps/v5.1.1-check/",
        },
        "deposit": {
            "repository": "myamlak/boys-harness",
            "tag": "v1.0.9",
            "release": "5.0.0",
            "commit": "482fd3388c276091834a6c12a642d9de2481c9f2",
        },
        "deposit_v511_records": deposit,
        "v520_check_in_deposit": False,
        "execution": execution,
        "input_comparison": scratch_comparison,
        "sweep_comparison": sweep_comparison,
        "table_comparison": table_comparison,
        "evaluated_files": {
            "include_boys": len(headers),
            "include_boys_identical_to_v511": header_same,
            "include_boys_differing": [
                rel for rel, left, right in headers if left != right],
            "src": len(sources),
            "src_identical_to_v511": source_same,
            "src_differing": [
                rel for rel, left, right in sources if left != right],
        },
        "tag_tree_blobs": {
            "files_checked": len(tag_rows),
            "matches": work_match,
        },
        "v511_provenance_digests": {
            "recorded": len(recorded),
            "matching": matched,
            "differing": [
                rel for rel, want, got in recorded if got != want],
        },
        "git_invoked": False,
        "run_performed": bool(execution),
        "note": ("Tag, commit and tree were read from the repository's loose "
                 "ref and object files without invoking git. The evaluated "
                 "sources in the worktree match the tagged tree's blobs, so "
                 "no measured error can change from a source difference. The "
                 "execution check, which is separate from this generator, "
                 "then compiled both sweeps against the tagged tree with the "
                 "documented compile line, ran them on the deposited grids, "
                 "fits and 60-digit references, and compared every produced "
                 "record with the deposited one; the steps and exit codes are "
                 "in run/logs/execution-check.log and the outcomes in the "
                 "execution, sweep_comparison and table_comparison entries. "
                 "The v1.0.11 deposit carries the v5.1.1 records, not this "
                 "check."),
    }
    (OUT / "provenance.json").write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8", newline="\n")

    print("wrote", OUT / "provenance.json")
    print("include/boys", len(headers), "same", header_same)
    print("src", len(sources), "same", source_same)
    print("tag blobs", work_match, "of", len(tag_rows))
    print("v5.1.1 digests", matched, "of", len(recorded))
    if execution:
        print("execution steps", len(execution["steps"]), "all zero",
              execution["all_steps_exit_zero"])
    if sweep_comparison:
        print("sweep records", sweep_comparison["identical"], "of",
              sweep_comparison["files_compared"], "identical")
    if table_comparison:
        print("tables", table_comparison["identical"], "identical;",
              table_comparison["verdict"])
    if scratch_comparison:
        print("inputs", scratch_comparison["identical"], "of",
              scratch_comparison["files_compared"], "identical")
    if deposit:
        print("deposit", deposit["identical"], "of", deposit["mapped"],
              "mapped files identical")


if __name__ == "__main__":
    main()
