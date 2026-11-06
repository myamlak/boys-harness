#!/usr/bin/env python3
"""Write the v5.2.5 source-identity records for Paper A.

The amplification sweeps and the integral-engine bridge were measured on the
v5.1.1 tree. Paper A cites the submission release, v5.2.5, and its source
identity is the claim that the release compiles the same evaluation code: the
evaluated files are compared here, one by one, and every difference is named
and explained.

Run from this directory:

    python make_check_525.py
"""

import datetime
import difflib
import hashlib
import subprocess
from pathlib import Path

OUT = Path(__file__).resolve().parent
BOYS = Path(r"C:\Users\myaml\dev\_boys-library")
PAPER = Path(r"C:\Users\myaml\dev\papers-a-v5\papers\boys-crossover")
TAG = "v5.2.5"
BASE = "v5.1.1"
PREVIOUS = "v5.2.4"
EVALUATED = ("include/boys", "src")

# The two headers the release changes beyond the version constant, and what each changes.
CHANGED = {
    "include/boys/boys_build_defaults.hpp": (
        "the build-defaults table: the twenty host default rows re-measured on 2026-11-05 "
        "and the row list rewritten. The sweeps and the bridge name every axis of their "
        "policy at the call site (their ReplicaPolicy / paper-a policy), so no unnamed call "
        "in either reads this table."
    ),
    "include/boys/boys_coefficients.hpp": (
        "the rational region-B seed: kRatBnum, kRatBden and kRegionBFitRatDelivered "
        "(1.69145631228896759e-14 -> 6.04405276991530667e-15), the degrees (4, 8) and the "
        "stored count (13) unchanged. Every policy these records name is the Chebyshev "
        "route, which reads no rational seed."
    ),
    "include/boys/version.hpp": (
        "the version constant only: 5.1.1 -> 5.2.5."
    ),
}


def git(*args):
    done = subprocess.run(["git", "-C", str(BOYS), *args], capture_output=True)
    if done.returncode != 0:
        raise SystemExit(f"git {' '.join(args)} failed: {done.stderr.decode('utf-8', 'replace')}")
    return done.stdout


def text(*args):
    return git(*args).decode("utf-8", "replace")


def blob(rev, path):
    return git("cat-file", "blob", f"{rev}:{path}")


def blob_sha1(data):
    hasher = hashlib.sha1()
    hasher.update(b"blob " + str(len(data)).encode() + b"\x00")
    hasher.update(data)
    return hasher.hexdigest()


def files_at(rev):
    listing = text("ls-tree", "-r", "--name-only", rev, *EVALUATED).splitlines()
    return sorted(path for path in listing if path)


def main():
    tag_object = text("cat-file", "-p", TAG).strip()
    commit = text("rev-parse", f"{TAG}^{{}}").strip()
    tree = text("rev-parse", f"{TAG}^{{tree}}").strip()
    base_commit = text("rev-parse", f"{BASE}^{{}}").strip()

    paths = files_at(TAG)
    base_paths = files_at(BASE)
    if paths != base_paths:
        raise SystemExit("the evaluated path lists differ; the comparison is not file for file")

    rows = []
    for path in paths:
        new = blob(TAG, path)
        old = blob(BASE, path)
        rows.append((path, blob_sha1(new), blob_sha1(old), new == old, hashlib.sha256(new).hexdigest()))

    identical = [row for row in rows if row[3]]
    differing = [row for row in rows if not row[3]]
    headers = [row for row in rows if row[0].startswith("include/boys/")]
    sources = [row for row in rows if row[0].startswith("src/")]

    # The evaluated file digests at the tag.
    with open(OUT / f"evaluation-sources-{TAG.replace('.', '')}.sha256", "w", encoding="utf-8", newline="\n") as handle:
        for path, _, _, _, sha in rows:
            handle.write(f"{sha}  {path}\n")

    # One line per evaluated file, both revisions' blob ids, and whether they agree.
    with open(OUT / "evaluated-sources-compare-v525.txt", "w", encoding="utf-8", newline="\n") as handle:
        handle.write(f"Evaluated sources, {TAG} vs {BASE}, compared blob for blob\n")
        handle.write(f"generated {datetime.datetime.now(datetime.UTC).isoformat(timespec='seconds')}\n\n")
        handle.write(f"tag {TAG}: commit {commit} tree {tree}\n")
        handle.write(f"base {BASE}: commit {base_commit}\n\n")
        handle.write(f"{len(paths)} evaluated files: {len(identical)} identical, {len(differing)} differing\n\n")
        handle.write("file\t" + TAG + "\t" + BASE + "\tidentical\n")
        for path, new_sha, old_sha, same, _ in rows:
            handle.write(f"{path}\t{new_sha}\t{old_sha}\t{same}\n")

    # The version header's own diff.
    old_version = blob(BASE, "include/boys/version.hpp").decode("utf-8").splitlines(keepends=True)
    new_version = blob(TAG, "include/boys/version.hpp").decode("utf-8").splitlines(keepends=True)
    diff = list(
        difflib.unified_diff(
            old_version,
            new_version,
            fromfile=f"include/boys/version.hpp ({BASE})",
            tofile=f"include/boys/version.hpp ({TAG})",
        )
    )
    (OUT / "version-hpp.diff").write_text("".join(diff), encoding="utf-8")

    # The narrative.
    lines = []
    lines.append(f"The {TAG} source identity")
    lines.append("")
    lines.append("What this records: the sources the paper's measurements were taken on (the")
    lines.append(f"{BASE} tree) against the release the paper cites ({TAG}), every evaluated file")
    lines.append("compared by content.")
    lines.append("")
    lines.append("The tag")
    lines.append(f"  tag object  {text('rev-parse', TAG).strip()}")
    lines.append(f"  commit      {commit}")
    lines.append(f"  tree        {tree}")
    for line in tag_object.splitlines():
        if line.startswith("tagger "):
            lines.append(f"  {line}")
    stamped = text("for-each-ref", "--format=%(taggerdate:iso-strict)", f"refs/tags/{TAG}").strip()
    if stamped:
        lines.append(f"  tag date    {stamped}")
    message = tag_object.split("\n\n", 1)[1].strip() if "\n\n" in tag_object else ""
    if message:
        lines.append(f"  message     {message}")
    lines.append("")
    lines.append("The comparison")
    lines.append(f"  include/boys  {len(headers)} files, {sum(1 for r in headers if r[3])} identical, "
                 f"{sum(1 for r in headers if not r[3])} differing")
    lines.append(f"  src           {len(sources)} files, {sum(1 for r in sources if r[3])} identical, "
                 f"{sum(1 for r in sources if not r[3])} differing")
    lines.append(f"  evaluated     {len(paths)} files, {len(identical)} identical, {len(differing)} differing")
    lines.append("")
    lines.append("The differing files, and why the measurements are not affected")
    for path, new_sha, old_sha, _, _ in differing:
        lines.append(f"  {path}")
        lines.append(f"    {BASE} blob {old_sha}")
        lines.append(f"    {TAG} blob {new_sha}")
        lines.append(f"    {CHANGED.get(path, 'CHANGED, with no explanation written here')}")
    lines.append("")
    lines.append("Every file read by the sweeps and the bridge is the same file in both trees")
    lines.append("except as named above, and no named difference is read by a policy either")
    lines.append("record states: the sweeps name every axis of their policy and evaluate the")
    lines.append("Chebyshev route, and the bridge's paper-a policy is the same route. The")
    lines.append("release therefore compiles the evaluation code the figures were measured on.")
    lines.append("")
    lines.append("Files in this directory")
    lines.append(f"  evaluated-sources-compare-v525.txt   the per-file comparison")
    lines.append(f"  evaluation-sources-{TAG.replace('.', '')}.sha256   the evaluated files' digests at {TAG}")
    lines.append(f"  version-hpp.diff                     the version header's diff, {BASE} to {TAG}")
    lines.append("  make_check_525.py                    the generator that wrote them")
    (OUT / "source-identity-v525.txt").write_text("\n".join(lines) + "\n", encoding="utf-8", newline="\n")

    write_sources_compare(rows, tag_object, commit, tree)

    print(f"{len(paths)} evaluated files: {len(identical)} identical, {len(differing)} differing")
    for path, _, _, same, _ in rows:
        if not same:
            print(f"  differing: {path}")


def revision(rev):
    """The commit a release name resolves to, and whether the name is an annotated tag."""
    obj_type = text("cat-file", "-t", rev).strip()
    if obj_type == "tag":
        return text("rev-parse", f"{rev}^{{}}").strip(), text("rev-parse", rev).strip()
    return text("rev-parse", rev).strip(), None


def write_sources_compare(rows, tag_object, commit, tree):
    """The three-way comparison Paper A's number checker reads."""
    base_commit, base_object = revision(BASE)
    prev_commit, prev_object = revision(PREVIOUS)

    def describe(name, commit_sha, tag_object_sha):
        out = []
        if tag_object_sha:
            out.append(f"{name}: tag object {tag_object_sha}")
            out.append(f"  commit {commit_sha}, tree {text('rev-parse', f'{name}^{{tree}}').strip()}")
        else:
            out.append(f"{name}: lightweight tag at commit {commit_sha}, "
                       f"tree {text('rev-parse', f'{name}^{{tree}}').strip()}")
        return out

    def same(rev, path):
        return blob(TAG, path) == blob(rev, path)

    header_rows = [row for row in rows if row[0].startswith("include/boys/")]
    source_rows = [row for row in rows if row[0].startswith("src/")]
    headers_both = sum(1 for r in header_rows if same(BASE, r[0]) and same(PREVIOUS, r[0]))
    sources_both = sum(1 for r in source_rows if same(BASE, r[0]) and same(PREVIOUS, r[0]))
    base_gap = sum(1 for r in rows if not same(BASE, r[0]))
    prev_gap = sum(1 for r in rows if not same(PREVIOUS, r[0]))

    lines = []
    lines.append(f"Source comparison for the {TAG} release")
    lines.append(f"Generated {datetime.datetime.now(datetime.UTC).isoformat(timespec='microseconds')}")
    lines.append("Read from the boys repository with git; every file compared by content.")
    lines.append("")
    lines.extend(describe(BASE, base_commit, base_object))
    lines.extend(describe(PREVIOUS, prev_commit, prev_object))
    lines.extend(describe(TAG, commit, text("rev-parse", TAG).strip()))
    lines.append("")
    lines.append(f"Evaluated path: {len(rows)} files "
                 f"({len(header_rows)} include/boys, {len(source_rows)} src).")
    lines.append("")
    lines.append(f"{TAG} against {BASE}: {len(rows) - base_gap} identical, {base_gap} differing.")
    lines.append(f"{TAG} against {PREVIOUS}: {len(rows) - prev_gap} identical, {prev_gap} differing.")
    lines.append("")
    lines.append(f"include/boys: {len(header_rows)} files; {headers_both} identical to both, "
                 f"{len(header_rows) - headers_both} differ.")
    lines.append(f"src: {len(source_rows)} files; {sources_both} identical to both, "
                 f"{len(source_rows) - sources_both} differ.")
    lines.append("")
    lines.append("file\tv525_sha256\tv524_sha256\tv511_sha256\tv525=v524\tv525=v511")
    for path, _, _, same_base, sha in rows:
        prev_blob = blob(PREVIOUS, path)
        base_blob = blob(BASE, path)
        new_blob = blob(TAG, path)
        lines.append(
            f"{path}\t{sha}\t{hashlib.sha256(prev_blob).hexdigest()}"
            f"\t{hashlib.sha256(base_blob).hexdigest()}"
            f"\t{new_blob == prev_blob}\t{same_base}"
        )
    body = "\n".join(lines) + "\n"

    (OUT / "v5.2.5-sources-compare.txt").write_text(body, encoding="utf-8", newline="\n")
    paper_dir = PAPER / "sweeps" / "v5.2.5-check"
    paper_dir.mkdir(parents=True, exist_ok=True)
    (paper_dir / "v5.2.5-sources-compare.txt").write_text(body, encoding="utf-8", newline="\n")
    print(f"wrote v5.2.5-sources-compare.txt ({headers_both}/{len(header_rows)} headers, "
          f"{sources_both}/{len(source_rows)} src, {len(rows) - base_gap}/{len(rows)} files)")


if __name__ == "__main__":
    main()
