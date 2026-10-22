"""Pull the 'Accuracy' section out of each rendered Boost.Math page.

BoostBook pages put every section in <h1..h6> headings whose text is the section title, so a
section runs from the end of its heading to the start of the next heading at the same or a
higher level.  The extracted text keeps table cells separated by ' | ' so the generated
accuracy tables stay readable as text.
"""

import html
import re
import sys
from pathlib import Path

HEADING = re.compile(r"<h([1-6])[^>]*>(.*?)</h\1>", re.S | re.I)


def text_of(fragment: str) -> str:
    s = re.sub(r"(?is)<(script|style).*?</\1>", " ", fragment)
    s = re.sub(r"(?i)<br\s*/?>", "\n", s)
    s = re.sub(r"(?i)</(p|div|li|h[1-6]|table|tr)>", "\n", s)
    s = re.sub(r"(?i)</t[dh]>", " | ", s)
    s = re.sub(r"<[^>]+>", "", s)
    s = html.unescape(s)
    s = re.sub(r"[ \t\r\f\v]+", " ", s)
    s = "\n".join(line.strip() for line in s.split("\n"))
    s = re.sub(r"\n{2,}", "\n", s)
    return s.strip()


def accuracy_sections(page: str):
    headings = [(m.start(), m.end(), int(m.group(1)), text_of(m.group(2)).strip())
                for m in HEADING.finditer(page)]
    out = []
    for i, (_, end, level, title) in enumerate(headings):
        if "accuracy" not in title.lower():
            continue
        stop = len(page)
        for (_, start, lvl, _) in headings[i + 1:]:
            if lvl <= level:
                stop = start
                break
        out.append((title, text_of(page[end:stop])))
    return out


def main() -> int:
    root = Path(sys.argv[1])
    outdir = Path(sys.argv[2])
    outdir.mkdir(parents=True, exist_ok=True)
    rows = []
    for path in sorted(root.rglob("*.html")):
        rel = path.relative_to(root).as_posix()
        page = path.read_text(encoding="utf-8", errors="replace")
        for n, (title, body) in enumerate(accuracy_sections(page)):
            dest = outdir / (rel.replace("/", "__") + f".acc{n}.txt")
            dest.write_text(body, encoding="utf-8")
            flags = []
            for token, label in ((r"\bulp", "ulp"), (r"epsilon", "epsilon"),
                                 (r"&epsilon;|ε|ϵ", "eps-glyph"),
                                 (r"\bMax\b", "Max"), (r"\bMean\b", "Mean")):
                if re.search(token, body, re.I | re.S):
                    flags.append(label)
            rows.append((rel, n, len(body), body.count("\n") + 1, ",".join(flags)))
    for rel, n, chars, lines, flags in rows:
        print(f"{rel}\tsec{n}\t{chars} chars\t{lines} lines\t{flags}")
    print(f"TOTAL pages with an Accuracy section: {len(rows)}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
