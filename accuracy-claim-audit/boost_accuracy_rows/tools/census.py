"""Count what each scoped accuracy section documents, in the shapes it documents it.

A section is a grid of cells: one per (test-data group, platform and floating-point type).
The cell's own number is the first "Max = ..." in it; the parenthesised ones are the same
measurement for other libraries, which are not claims about this library.  Cells whose content
begins with "(" report only other libraries and carry no claim of ours.  Everything the section
says in prose is reported separately, since prose is where a bound with a domain would live.
"""

import re
import sys
from pathlib import Path

MAXV = re.compile(r"Max\s*=\s*([+\-]?[\d.]+(?:e[+-]?\d+)?)\s*(?:&#949;|ε|ε)")
GROUP = re.compile(r"^[A-Z][^|]{3,70}$")


def cells(text: str):
    """Split an extracted section into cells, each a chunk between the '|' separators."""
    out, cur = [], []
    for line in text.split("\n"):
        if line.strip() == "|":
            out.append("\n".join(cur))
            cur = []
        else:
            cur.append(line)
    out.append("\n".join(cur))
    return [c for c in out if c.strip()]


def main() -> int:
    root = Path(sys.argv[1])
    for name in sys.argv[2:]:
        path = root / name
        text = path.read_text(encoding="utf-8")
        own, foreign, groups = [], 0, []
        for c in cells(text):
            m = MAXV.search(c)
            if not m:
                continue
            if c.strip().startswith("("):
                foreign += 1
                continue
            own.append(float(m.group(1)))
            # a group label line: the cell that holds only a name, before the number cells
            first = c.strip().split("\n")[0]
            if not re.match(r"^[+\-]?[\d.(]", first) and "=" not in first:
                groups.append(first)
        # group labels are separate cells in the grid; take them from the whole text
        for line in text.split("\n"):
            s = line.strip()
            if GROUP.match(s) and "=" not in s and not s.startswith(
                    ("GNU", "Sun", "Microsoft", "linux", "Win32", "Solaris", "Table", "Testing",
                     "Max", "Mean")):
                groups.append(s)
        print(f"{name}")
        print(f"   own max cells        : {len(own)}")
        print(f"   other-library cells  : {foreign}")
        print(f"   distinct group names : {len(set(groups))}")
        print(f"   own maxima           : {', '.join(f'{v:g}' for v in own)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
