#!/usr/bin/env python3
"""Compare the boys sources the bridge compiled with the v5.2.0 and v5.1.1 trees.

The bridge consumes external/boys as a plain copy, so what it compiled is
established by hashing that copy against the tagged worktrees.  No git command
is run.  The v5.2.0 worktree is the one verify_pub_tree.py checked; its path is
read from that script, its tag, commit and tree from pub-tree-identity.txt.

The comparison is against v5.2.0 (each file must be equal) and against v5.1.1
(where it shows which files the version bump changed).  The probe binaries are
dated against the newest file under external/boys, and the probe's own records
state the library version each run saw.
"""

import hashlib
import re
import sys
from datetime import datetime, timezone
from pathlib import Path

HERE = Path(__file__).resolve().parent
BRIDGE = Path(r"C:\Users\myaml\dev\bladawiec-bridge")
EXT = BRIDGE / "external" / "boys"
OLD = Path(r"C:\Users\myaml\dev\_boys-v511")
BUILDS = (BRIDGE / "build-bridge-default", BRIDGE / "build-bridge-paper-a")
RUNS = Path(r"C:\Users\myaml\dev\_bridge-v520\runs")
ARMS = ("default", "paper-a", "631gstar-default", "631gstar-paper-a",
        "ccpvtz-default", "ccpvtz-paper-a")


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def rel_files(root):
    if not root.is_dir():
        raise SystemExit("missing tree: %s" % root)
    return {p.relative_to(root).as_posix(): p
            for p in root.rglob("*") if p.is_file()}


def iso(path):
    return datetime.fromtimestamp(
        path.stat().st_mtime, timezone.utc).isoformat(timespec="seconds")


def main():
    out = Path(sys.argv[1]) if len(sys.argv) > 1 else None

    identity = (HERE / "pub-tree-identity.txt").read_text(encoding="utf-8")
    tag = re.search(r"refs/tags/(v[0-9.]+)", identity).group(1)
    commit = re.search(r"tag target commit\s+([0-9a-f]{40})", identity).group(1)
    tree = re.search(r"commit tree\s+([0-9a-f]{40})", identity).group(1)
    verifier = (HERE / "verify_pub_tree.py").read_text(encoding="utf-8")
    new_root = Path(re.search(r'^WORKTREE = Path\(r"(.*)"\)',
                              verifier, re.M).group(1))

    new = rel_files(new_root)
    old = rel_files(OLD)
    ext = rel_files(EXT)

    lines = [
        "Compiled-input comparison for the bridge tree, "
        + datetime.now().astimezone().isoformat(timespec="seconds"),
        "",
        "The tag's tree, read from the repository's object store:",
        "  tag     %s" % tag,
        "  commit  %s" % commit,
        "  tree    %s" % tree,
        "  worktree %s" % new_root,
        "",
        "bridge copy: %s" % EXT,
        "",
    ]

    checks = []

    def section(ext_paths):
        header_equal = old_equal = 0
        old_diff = []
        new_diff = []
        for rel in sorted(ext_paths):
            digest = sha256(EXT / rel)
            same_new = rel in new and sha256(new[rel]) == digest
            same_old = rel in old and sha256(old[rel]) == digest
            header_equal += same_new
            old_equal += same_old
            if not same_new:
                new_diff.append(rel)
            if not same_old:
                old_diff.append(rel)
        return header_equal, old_equal, new_diff, old_diff

    include_ext = [p for p in ext if p.startswith("include/boys/")]
    eq_new, eq_old, diff_new, diff_old = section(include_ext)
    lines.append(
        "include/boys: %d files; %d/%d equal %s and %d/%d equal v5.1.1."
        % (len(include_ext), eq_new, len(include_ext), tag,
           eq_old, len(include_ext)))
    if diff_new:
        lines.append("  differing from %s: %s" % (tag, ", ".join(diff_new)))
    if diff_old:
        lines.append("  differing from v5.1.1: %s" % ", ".join(diff_old))
    checks.append(eq_new == len(include_ext))

    src_ext = sorted(p for p in ext if p.startswith("src/"))
    lines.append("")
    lines.append("src/ files (all %d), bridge hash vs each worktree:" % len(src_ext))
    src_ok = True
    for rel in src_ext:
        digest = sha256(EXT / rel)
        same_new = rel in new and sha256(new[rel]) == digest
        same_old = rel in old and sha256(old[rel]) == digest
        src_ok = src_ok and same_new
        lines.append("%s\t%s\t%s %s %s %s" % (
            rel[len("src/"):], digest, tag, same_new, "v5.1.1", same_old))
    checks.append(src_ok)

    cmake = (BRIDGE / "integrals" / "CMakeLists.txt").read_text(encoding="utf-8")
    tused = []
    for name in re.findall(r"\$\{_qcxBoysExt\}/src/([A-Za-z0-9_.]+)", cmake):
        if name not in tused and not name.startswith("boys_cuda"):
            tused.append(name)
    lines.append("")
    lines.append(
        "compiled by qcx-integrals (integrals/CMakeLists.txt, CPU lane; the "
        "CUDA sources are behind the off CUDA lane):")
    lines.append("  " + ", ".join(tused))
    missing = [n for n in tused if "src/" + n not in src_ext]
    if missing:
        lines.append("  NOT under external/boys/src: " + ", ".join(missing))
    checks.append(len(tused) == 5 and not missing)

    newest = max(ext.values(), key=lambda p: p.stat().st_mtime)
    lines.append("")
    lines.append("builds and probe binaries:")
    exe_ok = True
    for build in BUILDS:
        cache = (build / "CMakeCache.txt").read_text(encoding="utf-8", errors="replace")
        generator = re.search(r"^CMAKE_GENERATOR:INTERNAL=(.*)$", cache, re.M)
        home = re.search(r"^CMAKE_HOME_DIRECTORY:INTERNAL=(.*)$", cache, re.M)
        build_type = re.search(r"^CMAKE_BUILD_TYPE:INTERNAL=(.*)$", cache, re.M)
        exe = build / "bridge" / "Release" / "qcx-bridge-probe.exe"
        lines.append("  %s" % build.name)
        lines.append("    generator  %s" % (generator.group(1) if generator else "?"))
        lines.append("    source     %s" % (home.group(1) if home else "?"))
        lines.append("    build type %s" % (
            build_type.group(1) if build_type else
            "multi-config, no CMAKE_BUILD_TYPE; Release used"))
        if not exe.is_file():
            lines.append("    MISSING %s" % exe)
            exe_ok = False
            continue
        newer = exe.stat().st_mtime > newest.stat().st_mtime
        exe_ok = exe_ok and newer
        lines.append("    binary     %s" % exe)
        lines.append("               %d bytes sha256 %s mtime %s"
                     % (exe.stat().st_size, sha256(exe), iso(exe)))
        lines.append("               newer than the newest external/boys file  %s"
                     % newer)
    lines.append("  newest external/boys file: %s  %s  %s" % (
        newest.relative_to(EXT).as_posix(), iso(newest), sha256(newest)))
    checks.append(exe_ok)

    lines.append("")
    lines.append("library version stated by each probe record (boys::VersionString):")
    version_ok = True
    for arm in ARMS:
        record = RUNS / (arm + ".boys_record.txt")
        if not record.is_file():
            lines.append("  %-20s MISSING" % arm)
            version_ok = False
            continue
        match = re.search(r"^\s+version\s+(\S+)\s+\(boys::VersionString\)",
                          record.read_text(encoding="utf-8", errors="replace"),
                          re.M)
        stated = match.group(1) if match else "?"
        version_ok = version_ok and stated == tag.lstrip("v")
        lines.append("  %-20s %s" % (arm, stated))
    checks.append(version_ok)

    lines.append("")
    lines.append("verdict: "
                 + ("bridge copy, binaries and records are the %s build" % tag
                    if all(checks) else "MISMATCH - see above"))
    text = "\n".join(lines) + "\n"
    if out is None:
        sys.stdout.write(text)
    else:
        out.write_text(text, encoding="utf-8", newline="\n")
    return 0 if all(checks) else 1


if __name__ == "__main__":
    sys.exit(main())
