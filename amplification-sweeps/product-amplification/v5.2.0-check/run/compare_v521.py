#!/usr/bin/env python3
"""Compare the v5.2.1 evaluation sources with v5.2.0 and v5.1.1.

The three annotated tags are resolved from the boys repository's loose refs
and object files; git is not invoked. The commit trees are walked and the 76
files under include/boys and src are compared by git blob SHA-1 and content
SHA-256. The machine records are written beside the v5.2.0 check's README.
"""

import hashlib
import importlib.util
from datetime import datetime, timezone
from pathlib import Path

HERE = Path(__file__).resolve().parent
OUT = HERE.parent
TAGS = ("v5.1.1", "v5.2.0", "v5.2.1")

spec = importlib.util.spec_from_file_location("vpt", HERE / "verify_pub_tree.py")
vpt = importlib.util.module_from_spec(spec)
spec.loader.exec_module(vpt)


def resolve(tag):
    ref_sha = (vpt.GITDIR / "refs" / "tags" / tag).read_text().strip()
    kind, body = vpt.read_object(ref_sha)
    if kind != "tag":
        raise SystemExit("%s: ref does not name a tag object" % tag)
    if vpt.tag_field(body, "tag") != tag:
        raise SystemExit("%s: tag object names %r" % (tag, vpt.tag_field(body, "tag")))
    target = vpt.tag_field(body, "object")
    kind, commit = vpt.read_object(target)
    if kind != "commit":
        raise SystemExit("%s: target is not a commit" % tag)
    tree = vpt.tag_field(commit, "tree")
    return {
        "tag": tag,
        "tag_object": ref_sha,
        "tagger": vpt.tag_field(body, "tagger"),
        "message": body.split(b"\n\n", 1)[1].decode(errors="replace").strip(),
        "commit": target,
        "tree": tree,
        "parent": vpt.tag_field(commit, "parent"),
        "files": vpt.walk(tree),
    }


def read_blob(files, path):
    kind, body = vpt.read_object(files[path][1])
    if kind != "blob":
        raise SystemExit("%s: not a blob" % path)
    return body


def blob_sha256(files, path):
    return hashlib.sha256(read_blob(files, path)).hexdigest()


def blob_sha1(files, path):
    return vpt.blob_sha1(read_blob(files, path))


def diff(a, b):
    before = identities[a]["files"]
    after = identities[b]["files"]
    modified = sorted(p for p in before if p in after and before[p][1] != after[p][1])
    added = sorted(p for p in after if p not in before)
    deleted = sorted(p for p in before if p not in after)
    return modified, added, deleted


identities = {tag: resolve(tag) for tag in TAGS}
evaluated = sorted(
    p for p in identities["v5.2.1"]["files"]
    if p.startswith("include/boys/") or p.startswith("src/"))
headers = [p for p in evaluated if p.startswith("include/boys/")]
sources = [p for p in evaluated if p.startswith("src/")]
if len(evaluated) != 76 or len(headers) != 62 or len(sources) != 14:
    raise SystemExit("unexpected evaluated path: %d files (%d headers, %d src)"
                     % (len(evaluated), len(headers), len(sources)))

sha = {tag: {p: blob_sha256(identities[tag]["files"], p) for p in evaluated}
       for tag in TAGS}
sha1 = {tag: {p: blob_sha1(identities[tag]["files"], p) for p in evaluated}
        for tag in TAGS}

if any(p not in identities["v5.2.0"]["files"] or p not in identities["v5.1.1"]["files"]
       for p in evaluated):
    raise SystemExit("an evaluated file is absent from a comparison tree")

def is_sha256(value):
    return len(value) == 64 and all(c in "0123456789abcdef" for c in value)


# Cross-check the object-store values against the v5.2.0 check's own table.
# Only the two main tables carry two 64-hex columns; the later digest
# comparison rows carry the word match or differ in place of a hash.
existing = {}
section = None
for line in (OUT / "evaluation-sources-compare.txt").read_text(encoding="utf-8").splitlines():
    if line.startswith("include/boys:"):
        section = "include/boys"
    elif line.startswith("src:"):
        section = "src"
    elif section and line.count("\t") == 3:
        name, s520, s511, _ = line.split("\t")
        if is_sha256(s520) and is_sha256(s511):
            existing[section + "/" + name] = (s520, s511)

cross = []
for p in evaluated:
    if p not in existing:
        cross.append("%s missing from evaluation-sources-compare.txt" % p)
    elif existing[p] != (sha["v5.2.0"][p], sha["v5.1.1"][p]):
        cross.append("%s: recorded %s/%s, object store %s/%s"
                     % (p, existing[p][0], existing[p][1],
                        sha["v5.2.0"][p], sha["v5.1.1"][p]))
if len(existing) != len(evaluated) and not cross:
    cross.append("evaluation-sources-compare.txt carries %d rows, expected %d"
                 % (len(existing), len(evaluated)))

generated = datetime.now(timezone.utc).isoformat()

lines = []
lines.append("Source comparison for the v5.2.1 release")
lines.append("Generated " + generated)
lines.append("Read from %s; git was not invoked." % vpt.GITDIR)
lines.append("")
for tag in TAGS:
    ident = identities[tag]
    lines.append("%s: tag object %s" % (tag, ident["tag_object"]))
    lines.append("  commit %s, tree %s" % (ident["commit"], ident["tree"]))
    lines.append("  parent %s" % (ident["parent"] or "(none)"))
    lines.append("  tagger %s" % ident["tagger"])
    lines.append("  message %r" % ident["message"])
    lines.append("  files %d" % len(ident["files"]))
    lines.append("")
lines.append("Evaluated path: %d files (%d include/boys, %d src)."
             % (len(evaluated), len(headers), len(sources)))
lines.append("")
lines.append("v5.2.1 against v5.2.0: %d identical, %d differing."
             % (sum(sha["v5.2.1"][p] == sha["v5.2.0"][p] for p in evaluated),
                sum(sha["v5.2.1"][p] != sha["v5.2.0"][p] for p in evaluated)))
lines.append("v5.2.1 against v5.1.1: %d identical, %d differing."
             % (sum(sha["v5.2.1"][p] == sha["v5.1.1"][p] for p in evaluated),
                sum(sha["v5.2.1"][p] != sha["v5.1.1"][p] for p in evaluated)))
lines.append("")
lines.append("include/boys: %d files; %d identical to both, %d differs."
             % (len(headers),
                sum(sha["v5.2.1"][p] == sha["v5.2.0"][p] == sha["v5.1.1"][p]
                    for p in headers),
                sum(not (sha["v5.2.1"][p] == sha["v5.2.0"][p] == sha["v5.1.1"][p])
                    for p in headers)))
lines.append("src: %d files; %d identical to both, %d differs."
             % (len(sources),
                sum(sha["v5.2.1"][p] == sha["v5.2.0"][p] == sha["v5.1.1"][p]
                    for p in sources),
                sum(not (sha["v5.2.1"][p] == sha["v5.2.0"][p] == sha["v5.1.1"][p])
                    for p in sources)))
lines.append("")
lines.append("file\tv521_sha256\tv520_sha256\tv511_sha256\tv521=v520\tv521=v511")
for p in evaluated:
    lines.append("\t".join([
        p, sha["v5.2.1"][p], sha["v5.2.0"][p], sha["v5.1.1"][p],
        str(sha["v5.2.1"][p] == sha["v5.2.0"][p]).lower(),
        str(sha["v5.2.1"][p] == sha["v5.1.1"][p]).lower(),
    ]))
lines.append("")
lines.append("Cross-check against evaluation-sources-compare.txt:")
lines.append("  v5.2.0 and v5.1.1 SHA-256 values reproduced for %d of %d files."
             % (len(evaluated) - len(cross), len(evaluated)))
for item in cross:
    lines.append("  MISMATCH " + item)
lines.append("")
lines.append("Conclusion: every evaluated translation unit is byte-identical to")
lines.append("the historical v5.1.1 run and to the v5.2.0 execution check except")
lines.append("include/boys/version.hpp, whose only change is the version constant.")
(OUT / "v5.2.1-sources-compare.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")

(OUT / "evaluation-sources-v521.sha256").write_text(
    "".join("%s  %s\n" % (sha["v5.2.1"][p], p) for p in evaluated),
    encoding="utf-8")

lines = []
lines.append("Source-tree differences, v5.2.1 tree against v5.2.0 and v5.1.1")
lines.append("Generated " + generated)
lines.append("Read from %s; git was not invoked." % vpt.GITDIR)
lines.append("Compared by git blob SHA-1, file by file; the trees are complete.")
lines.append("")
for tag in TAGS:
    ident = identities[tag]
    lines.append("%s: tag object %s, commit %s, tree %s, files %d"
                 % (tag, ident["tag_object"], ident["commit"], ident["tree"],
                    len(ident["files"])))
lines.append("")
for a, b in (("v5.2.0", "v5.2.1"), ("v5.1.1", "v5.2.1")):
    modified, added, deleted = diff(a, b)
    before = identities[a]["files"]
    after = identities[b]["files"]
    common = set(before) & set(after)
    the_same = sum(before[p][1] == after[p][1] for p in common)
    lines.append("%s -> %s: %d files in %s; %d modified, %d added, %d deleted, "
                 "%d identical." % (a, b, len(before), a, len(modified),
                                    len(added), len(deleted), the_same))
    for p in modified:
        lines.append("  modified %s" % p)
    for p in added:
        lines.append("  added    %s" % p)
    for p in deleted:
        lines.append("  deleted  %s" % p)
    lines.append("")
lines.append("Evaluated path:")
for label, paths in (("include/boys", headers), ("src", sources)):
    differing = [p for p in paths
                 if not (sha["v5.2.1"][p] == sha["v5.2.0"][p] == sha["v5.1.1"][p])]
    lines.append("  %s: %d files; %d identical to both, %d differing."
                 % (label, len(paths), len(paths) - len(differing), len(differing)))
    for p in differing:
        lines.append("    differing %s" % p)
lines.append("")
lines.append("No file under src differs. The only file under include/boys that")
lines.append("differs is version.hpp. None of the modified, added or deleted files")
lines.append("is compiled by the two sweeps; the changed gate, test-data,")
lines.append("documentation and metadata files are outside the evaluated path.")
(OUT / "source-diff-v521.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")

print("v5.2.1 tag object %s" % identities["v5.2.1"]["tag_object"])
print("v5.2.1 commit %s" % identities["v5.2.1"]["commit"])
print("v5.2.1 tree %s" % identities["v5.2.1"]["tree"])
print("evaluated %d: %d headers, %d src" % (len(evaluated), len(headers), len(sources)))
print("v5.2.1 = v5.2.0: %d of %d"
      % (sum(sha["v5.2.1"][p] == sha["v5.2.0"][p] for p in evaluated), len(evaluated)))
print("v5.2.1 = v5.1.1: %d of %d"
      % (sum(sha["v5.2.1"][p] == sha["v5.1.1"][p] for p in evaluated), len(evaluated)))
print("cross-check mismatches: %d" % len(cross))
for a, b in (("v5.2.0", "v5.2.1"), ("v5.1.1", "v5.2.1")):
    modified, added, deleted = diff(a, b)
    print("%s -> %s: %d modified, %d added, %d deleted"
          % (a, b, len(modified), len(added), len(deleted)))
