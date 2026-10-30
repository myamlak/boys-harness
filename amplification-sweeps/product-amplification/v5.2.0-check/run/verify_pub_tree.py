#!/usr/bin/env python3
"""Verify the clean v5.2.0 worktree against the tag's own tree, without git.

The annotated tag object, its target commit and the commit's tree are read
from the repository's object store (loose files and packfiles).  Every file
the tree records is compared with the working tree by git blob SHA-1, and
every working-tree file is checked against the tree, so an untracked or
modified file is reported.

This stands in for `git describe --exact-match` and `git status --porcelain`:
it establishes the same two facts (HEAD is the tag's commit; the tree is
clean) from the files on disk.  Exits non-zero on any difference.
"""

import hashlib
import sys
import zlib
from functools import lru_cache
from pathlib import Path

GITDIR = Path(r"C:\Users\myaml\dev\_mainonly\.git")
WORKTREE = Path(r"C:\Users\myaml\dev\_boys-v520-pub")
WORKTREE_GITDIR = GITDIR / "worktrees" / "_boys-v520-pub"
TAG_NAME = "v5.2.0"

OUT = Path(__file__).resolve().parent / "pub-tree-identity.txt"

TYPES = {1: "commit", 2: "tree", 3: "blob", 4: "tag"}


class MissingObject(Exception):
    def __init__(self, sha):
        super().__init__(sha)
        self.sha = sha


@lru_cache(maxsize=None)
def pack_bytes(path):
    return Path(path).read_bytes()


@lru_cache(maxsize=None)
def pack_index(path):
    data = Path(path).read_bytes()
    if data[:4] != b"\xfftOc":
        raise SystemExit("unsupported idx: %s" % path)
    version = int.from_bytes(data[4:8], "big")
    if version != 2:
        raise SystemExit("unsupported idx version %d: %s" % (version, path))
    fanout = [int.from_bytes(data[8 + 4 * i:12 + 4 * i], "big") for i in range(256)]
    count = fanout[255]
    names = 8 + 1024
    offsets = names + 20 * count + 4 * count
    large = offsets + 4 * count
    return data, fanout, names, offsets, large, count


def idx_offset(idx_path, sha_hex):
    data, fanout, names, offsets, large, count = pack_index(idx_path)
    target = bytes.fromhex(sha_hex)
    bucket = target[0]
    lo = fanout[bucket - 1] if bucket else 0
    hi = fanout[bucket]
    while lo < hi:
        mid = (lo + hi) // 2
        name = data[names + 20 * mid:names + 20 * mid + 20]
        if name < target:
            lo = mid + 1
        elif name > target:
            hi = mid
        else:
            offset = int.from_bytes(data[offsets + 4 * mid:offsets + 4 * mid + 4], "big")
            if offset & 0x80000000:
                index = offset & 0x7FFFFFFF
                offset = int.from_bytes(data[large + 8 * index:large + 8 * index + 8], "big")
            return offset
    return None


def inflate(data, offset):
    engine = zlib.decompressobj()
    body = engine.decompress(data[offset:])
    consumed = len(data) - offset - len(engine.unused_data)
    return body, offset + consumed


def header(data, offset):
    byte = data[offset]
    offset += 1
    kind = (byte >> 4) & 7
    size = byte & 15
    shift = 4
    while byte & 0x80:
        byte = data[offset]
        offset += 1
        size |= (byte & 0x7F) << shift
        shift += 7
    return kind, size, offset


def varint(data, pos):
    value = 0
    shift = 0
    while True:
        byte = data[pos]
        pos += 1
        value |= (byte & 0x7F) << shift
        shift += 7
        if not byte & 0x80:
            return value, pos


def apply_delta(base, delta):
    base_size, pos = varint(delta, 0)
    result_size, pos = varint(delta, pos)
    if base_size != len(base):
        raise SystemExit("delta base size mismatch")
    out = bytearray()
    while pos < len(delta):
        op = delta[pos]
        pos += 1
        if op & 0x80:
            copy_offset = 0
            copy_size = 0
            for i in range(4):
                if op & (1 << i):
                    copy_offset |= delta[pos] << (8 * i)
                    pos += 1
            for i in range(3):
                if op & (1 << (4 + i)):
                    copy_size |= delta[pos] << (8 * i)
                    pos += 1
            if copy_size == 0:
                copy_size = 0x10000
            out += base[copy_offset:copy_offset + copy_size]
        elif op:
            out += delta[pos:pos + op]
            pos += op
        else:
            raise SystemExit("bad delta opcode")
    if len(out) != result_size:
        raise SystemExit("delta result size mismatch")
    return bytes(out)


def read_pack_at(path, data, offset, depth=0):
    if depth > 64:
        raise SystemExit("delta chain too deep")
    kind, size, pos = header(data, offset)
    if kind in TYPES:
        body, _ = inflate(data, pos)
        if len(body) != size:
            raise SystemExit("pack size mismatch")
        return TYPES[kind], body
    if kind == 6:
        byte = data[pos]
        pos += 1
        relative = byte & 0x7F
        while byte & 0x80:
            byte = data[pos]
            pos += 1
            relative = ((relative + 1) << 7) | (byte & 0x7F)
        base_kind, base = read_pack_at(path, data, offset - relative, depth + 1)
        delta, _ = inflate(data, pos)
        return base_kind, apply_delta(base, delta)
    if kind == 7:
        base_sha = data[pos:pos + 20].hex()
        pos += 20
        base_kind, base = read_object(base_sha)
        delta, _ = inflate(data, pos)
        return base_kind, apply_delta(base, delta)
    raise SystemExit("bad pack object type %d" % kind)


def read_object(sha):
    loose = GITDIR / "objects" / sha[:2] / sha[2:]
    if loose.is_file():
        raw = zlib.decompress(loose.read_bytes())
        head, _, body = raw.partition(b"\0")
        kind, _, size = head.partition(b" ")
        if int(size) != len(body):
            raise SystemExit("loose object size mismatch: %s" % sha)
        return kind.decode(), body
    for idx in sorted((GITDIR / "objects" / "pack").glob("*.idx")):
        pack = idx.with_suffix(".pack")
        if not pack.is_file():
            continue
        offset = idx_offset(idx, sha)
        if offset is None:
            continue
        return read_pack_at(pack, pack_bytes(pack), offset)
    raise MissingObject(sha)


def parse_tree(body):
    entries = {}
    offset = 0
    while offset < len(body):
        end = body.index(b"\0", offset)
        mode, name = body[offset:end].split(b" ", 1)
        entry = body[end + 1:end + 21].hex()
        entries[name.decode()] = (mode.decode(), entry)
        offset = end + 21
    return entries


def walk(sha, prefix=""):
    kind, body = read_object(sha)
    if kind != "tree":
        raise SystemExit("not a tree: %s (%s)" % (sha, kind))
    files = {}
    for name, (mode, entry) in parse_tree(body).items():
        path = prefix + name
        if mode == "40000":
            files.update(walk(entry, path + "/"))
        else:
            files[path] = (mode, entry)
    return files


def blob_sha1(data):
    return hashlib.sha1(b"blob %d\0" % len(data) + data).hexdigest()


def tag_field(body, name):
    prefix = name.encode() + b" "
    for line in body.split(b"\n"):
        if line.startswith(prefix):
            return line[len(prefix):].decode()
    return None


def main():
    tag_object = (GITDIR / "refs" / "tags" / TAG_NAME).read_text().strip()
    kind, tag_body = read_object(tag_object)
    if kind != "tag":
        raise SystemExit("tag ref does not name a tag object: %s" % kind)
    if tag_field(tag_body, "tag") != TAG_NAME:
        raise SystemExit("tag object names %r" % tag_field(tag_body, "tag"))
    target = tag_field(tag_body, "object")

    kind, commit_body = read_object(target)
    if kind != "commit":
        raise SystemExit("tag target is not a commit: %s" % kind)
    tree = tag_field(commit_body, "tree")
    if tree is None:
        raise SystemExit("commit carries no tree line")

    head = (WORKTREE_GITDIR / "HEAD").read_text().strip()

    tracked = walk(tree)
    missing = []
    differing = []
    checked = 0
    for relative, (mode, expected) in sorted(tracked.items()):
        path = WORKTREE / relative
        if not path.is_file():
            missing.append(relative)
            continue
        actual = blob_sha1(path.read_bytes())
        checked += 1
        if actual != expected:
            differing.append((relative, expected, actual))

    working = set()
    for path in WORKTREE.rglob("*"):
        if path.is_file():
            relative = path.relative_to(WORKTREE).as_posix()
            if relative != ".git":
                working.add(relative)
    extra = sorted(working - set(tracked))

    report = []
    report.append("v5.2.0 source tree identity, read from the repository's object store")
    report.append("  tag ref refs/tags/%s -> tag object %s" % (TAG_NAME, tag_object))
    report.append("  tag target commit          %s" % target)
    report.append("  commit tree                %s" % tree)
    report.append("  worktree HEAD              %s" % head)
    report.append("  HEAD is the tag's commit   %s" % (head == target))
    report.append("  tracked files in tree      %d" % len(tracked))
    report.append("  checked against worktree   %d" % checked)
    report.append("  differing                  %d" % len(differing))
    report.append("  missing from worktree      %d" % len(missing))
    report.append("  extra in worktree          %d" % len(extra))
    for item in differing:
        report.append("  DIFF %s expected %s actual %s" % item)
    for item in missing:
        report.append("  MISSING %s" % item)
    for item in extra:
        report.append("  EXTRA %s" % item)

    ok = head == target and not differing and not missing and not extra
    report.append("  verdict                    %s" % ("CLEAN AT THE TAG" if ok else "NOT CLEAN"))

    text = "\n".join(report) + "\n"
    OUT.write_text(text, encoding="utf-8")
    print(text, end="")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
