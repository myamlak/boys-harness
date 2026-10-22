"""Extract candidate accuracy statements from the crlibm manual.

Strips LaTeX comments and markup, splits into sentences, keeps those that
assert something about accuracy.  Classification is done by hand afterwards;
this script only produces the candidate list with file:line provenance.
"""
import os
import re
import sys

DOC = r"C:\Users\myaml\dev\_claim-audit-crlibm\src\crlibm\docs\latex"

VOCAB = re.compile(
    r"accuracy|correctly rounded|correct rounding|worst case|worst-case|"
    r"error bound|ulp|guarantee|faithful", re.I)

# sections that are pure toolkit / proof-mechanics discussion
SKIP_FILES = set()


def strip_tex(s):
    s = re.sub(r"(?<!\\)%.*$", "", s)            # comments
    s = re.sub(r"\\(label|ref|pageref|cite|citeauthor)\{[^}]*\}", " ", s)
    s = re.sub(r"\\[a-zA-Z]+\*?(\[[^\]]*\])?", " ", s)
    s = s.replace("$", " ").replace("{", " ").replace("}", " ")
    s = s.replace("\\", " ").replace("~", " ")
    s = re.sub(r"\s+", " ", s)
    return s.strip()


def sentences(path):
    """Yield (lineno, sentence) for the whole file."""
    text_lines = []
    with open(path, encoding="utf-8", errors="replace") as fh:
        for i, line in enumerate(fh, 1):
            text_lines.append((i, line))
    buf, start = "", None
    for i, line in text_lines:
        clean = strip_tex(line)
        if start is None and clean:
            start = i
        buf += " " + clean
        if clean.endswith(".") or clean == "":
            for s in re.split(r"(?<=[.:;])\s+", buf):
                s = s.strip()
                if s:
                    yield start, s
            buf, start = "", None
    if buf.strip():
        yield start or 0, buf.strip()


def main():
    out = []
    for fn in sorted(os.listdir(DOC)):
        if not fn.endswith(".tex") or fn in SKIP_FILES:
            continue
        for lineno, s in sentences(os.path.join(DOC, fn)):
            if len(s) < 15:
                continue
            if VOCAB.search(s):
                out.append((fn, lineno, s))
    with open("claim_candidates.txt", "w", encoding="utf-8") as fh:
        for fn, lineno, s in out:
            fh.write("%s:%d\t%s\n" % (fn, lineno, s))
    print("candidate accuracy sentences:", len(out))
    per = {}
    for fn, _, _ in out:
        per[fn] = per.get(fn, 0) + 1
    for fn in sorted(per, key=lambda k: -per[k]):
        print("  %-32s %d" % (fn, per[fn]))


if __name__ == "__main__":
    main()
