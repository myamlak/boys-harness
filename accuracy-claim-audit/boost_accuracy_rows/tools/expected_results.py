"""Parse add_expected_result(...) calls out of the Boost.Math test suite.

Each call is one asserted bound: a compiler / stdlib / platform / type / data-group
selector (all regexes) together with the maximum peak and mean error, in units of
epsilon, that the test suite will accept for the group.  These are the library's
machine-checked accuracy bounds; the tables in the reference manual are the maxima
observed on one run.
"""

import csv
import re
import sys
from pathlib import Path

CALL = re.compile(r"add_expected_result\s*\(", re.S)


def strip_comments(text: str) -> str:
    """Drop // to end of line, and /* */, outside string literals.

    The calls carry a comment explaining each argument, so an argument and the
    comment belonging to the one before it would otherwise run together.
    """
    out, i, in_str, in_block = [], 0, False, False
    while i < len(text):
        c = text[i]
        if in_str:
            out.append(c)
            if c == "\\":
                out.append(text[i + 1])
                i += 2
                continue
            if c == '"':
                in_str = False
        elif in_block:
            if c == "*" and text[i + 1:i + 2] == "/":
                in_block = False
                i += 2
                continue
        elif c == '"':
            in_str = True
            out.append(c)
        elif c == "/" and text[i + 1:i + 2] == "/":
            while i < len(text) and text[i] != "\n":
                i += 1
            continue
        elif c == "/" and text[i + 1:i + 2] == "*":
            in_block = True
            i += 2
            continue
        else:
            out.append(c)
        i += 1
    return "".join(out)


DEFINED = {"BOOST_MATH_ENABLE_SYCL": False}
DIRECTIVES = set()


def evaluate_directive(line: str, macros):
    """Decide a preprocessor conditional inside an argument list.

    Some calls carry their trailing arguments in #ifdef branches, so which bound a
    configuration gets depends on the macros.  Only the branch the configuration takes may be
    read; the others are dropped.  Every directive met is recorded, so the report can say which
    conditions the numbers stood behind.
    """
    body = line.lstrip("#").strip()
    name = body.split()[0] if body.split() else body
    macros.add(name + (": " + " ".join(body.split()[1:]) if len(body.split()) > 1 else ""))
    if name in ("ifdef", "ifndef", "if"):
        if name == "ifdef":
            return DEFINED.get(body.split()[1], False) if len(body.split()) > 1 else True
        if name == "ifndef":
            return not DEFINED.get(body.split()[1], False) if len(body.split()) > 1 else False
        return True
    return None  # else / endif


def split_args(text: str, macros=None):
    """Return the top-level comma separated arguments of a call whose '(' was consumed."""
    if macros is None:
        macros = DIRECTIVES
    args, depth, cur, i, in_str = [], 1, "", 0, False
    in_dropped = False  # a conditional branch whose arguments do not apply
    while i < len(text):
        c = text[i]
        if not in_str and c == "#":
            end = text.find("\n", i)
            line = text[i:end if end >= 0 else len(text)]
            verdict = evaluate_directive(line, macros)
            if verdict is None:
                if line.strip().startswith("#else"):
                    in_dropped = not in_dropped
                elif line.strip().startswith("#endif"):
                    in_dropped = False
            else:
                if verdict:
                    in_dropped = False
                else:
                    in_dropped = True
            i = (end + 1) if end >= 0 else len(text)
            continue
        if in_dropped:
            i += 1
            continue
        if in_str:
            cur += c
            if c == "\\":
                cur += text[i + 1]
                i += 2
                continue
            if c == '"':
                in_str = False
        elif c == '"':
            in_str = True
            cur += c
        elif c == "(":
            depth += 1
            cur += c
        elif c == ")":
            depth -= 1
            if depth == 0:
                args.append(cur)
                return args, i + 1
            cur += c
        elif c == "," and depth == 1:
            args.append(cur)
            cur = ""
        else:
            cur += c
        i += 1
    return args, i


def unquote(raw: str) -> str:
    raw = raw.strip()
    m = re.fullmatch(r'"(.*)"', raw, re.S)
    return m.group(1) if m else raw


def main() -> int:
    root = Path(sys.argv[1])
    out = Path(sys.argv[2])
    rows = []
    macros = set()
    for path in sorted(root.glob("test_*.cpp")):
        text = strip_comments(path.read_text(encoding="utf-8", errors="replace"))
        for m in CALL.finditer(text):
            args, _ = split_args(text[m.end():], macros)
            if len(args) < 8:
                continue
            line = text[:m.start()].count("\n") + 1
            rows.append([path.name, line] + [unquote(a).replace("\n", " ") for a in args[:9]])
    with out.open("w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["file", "line", "compiler", "stdlib", "platform", "type", "group",
                    "function", "max_eps", "mean_eps"])
        w.writerows(rows)
    print(f"files scanned : {len(list(root.glob('test_*.cpp')))}")
    print(f"bound records : {len(rows)}")
    print("directives met: " + ", ".join(sorted(macros)))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
