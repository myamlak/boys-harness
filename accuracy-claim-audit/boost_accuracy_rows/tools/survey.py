"""Two lookups the row table needs.

1. Domain of a named data group: the Boost.Math data files (test/*.ipp) hold arrays of
   published test points, and the group labels in the reference manual are the names of those
   arrays.  This prints the range of each column so a claim's domain can be stated as an
   interval with something behind it.

2. The asserted bound for one configuration: the test suite selects the first
   add_expected_result whose six regexes match the string
   compiler|stdlib|platform|type|data-group|function, and falls back to (1, 1) ulp when none
   matches.  This finds that bound for a configuration given on the command line.
"""

import csv
import re
import sys
from pathlib import Path

SC = re.compile(r"SC_\(([^)]*)\)")
ARRAY = re.compile(r"std::array<std::array<[^,]+,\s*(\d+)>,\s*(\d+)>\s+(\w+)\s*=\s*\{\s*\{")
NUMBER = re.compile(r"^[+-]?[\d.]+(?:e[+-]?\d+)?$")


def arrays_in(path: Path):
    text = path.read_text(encoding="utf-8", errors="replace")
    for m in ARRAY.finditer(text):
        name = m.group(3)
        cols = int(m.group(1))
        rows = int(m.group(2))
        body = text[m.end():]
        values, i = [], 0
        for mm in SC.finditer(body):
            raw = mm.group(1).strip()
            if not NUMBER.match(raw):
                break
            values.append(float(raw))
            i += 1
            if i == cols * rows:
                break
        if len(values) < cols * rows:
            continue
        grid = [values[j * cols:(j + 1) * cols] for j in range(rows)]
        yield name, cols, grid


def report_domains(test_dir: Path, names):
    for path in sorted(test_dir.glob("*.ipp")):
        for name, cols, grid in arrays_in(path):
            if names and name not in names:
                continue
            lo = [min(r[c] for r in grid) for c in range(cols)]
            hi = [max(r[c] for r in grid) for c in range(cols)]
            span = "  ".join(f"c{c}: [{lo[c]:.6g}, {hi[c]:.6g}]" for c in range(cols))
            print(f"{path.name:28s} {name:22s} n={len(grid):4d}  {span}")


LARGEST_TYPE = r"(long\s+)?double|real_concept"


def matches(record, config):
    """Replicate the suite's selection: every selector regex must match the whole field.

    The type selector is not always a literal: the tests pass a C++ variable named
    largest_type, whose value is a regex covering the widest floating-point type this
    platform has.  Where long double is no wider than double, that value is
    (long\\s+)?double|real_concept, which is what a double build matches here.
    """
    fields = ("compiler", "stdlib", "platform", "type", "group", "function")
    for k in fields:
        pattern = record[k]
        if k == "type" and pattern.strip() == "largest_type":
            pattern = LARGEST_TYPE
        try:
            if not re.fullmatch(pattern, config[k], re.I):
                return False
        except re.error:
            return False
    return True


def report_bound(bounds_csv: Path, file_name, config):
    rows = [r for r in csv.DictReader(bounds_csv.open(encoding="utf-8"))
            if r["file"] == file_name]
    for r in rows:
        if matches(r, config):
            print(f"first match in {file_name}: line {r['line']}  "
                  f"group='{r['group']}' function='{r['function']}' "
                  f"type='{r['type']}' platform='{r['platform']}' -> {r['max_eps']} eps")
            return
    print(f"no match in {file_name}: the suite would use its default of (1, 1) eps")


if __name__ == "__main__":
    if sys.argv[1] == "domains":
        report_domains(Path(sys.argv[2]), set(sys.argv[3:]))
    else:
        config = dict(zip(("compiler", "stdlib", "platform", "type", "group", "function"),
                          sys.argv[3:9]))
        report_bound(Path(sys.argv[2]), sys.argv[1], config)
