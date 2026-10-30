#!/usr/bin/env python3
"""Collect the bridge lane's FINISH lines and commands into one record.

The lane script appends one FINISH line per step to logs/<step>.log.  This
script reads those lines, pairs them with the commands exactly as the lane
script writes them, and resolves the path variables the commands use.  A step
is reported as it finished; no value is restated from memory.
"""

import re
import sys
from datetime import datetime
from pathlib import Path

LOGS = Path(r"C:\Users\myaml\dev\_bridge-v520\logs")
SCRIPT = Path(r"C:\Users\myaml\dev\_bridge-v520\scripts\run_all.ps1")

FINISH = re.compile(r"^FINISH\s+(\S+)\s+STEP=(\S+)\s+EXIT=(\d+)\s+SECONDS=(\S+)$")
ASSIGN = re.compile(r'^\$(\w+)\s*=\s*[\'"](.*)[\'"]$')


def main():
    out = Path(sys.argv[1]) if len(sys.argv) > 1 else None

    steps = []
    for log in sorted(LOGS.glob("*.log")):
        for line in log.read_text(encoding="utf-8", errors="replace").splitlines():
            match = FINISH.match(line.strip())
            if match:
                steps.append((match.group(1), match.group(2), int(match.group(3)),
                              match.group(4), log.name, line.strip()))
    steps.sort()

    text = SCRIPT.read_text(encoding="utf-8", errors="replace")
    values = []
    commands = []
    for line in text.splitlines():
        stripped = line.strip()
        assign = ASSIGN.match(stripped)
        if assign:
            values.append((assign.group(1), assign.group(2)))
        elif stripped.startswith("& "):
            commands.append(stripped)

    # Only the variables the commands reach, directly or through another
    # value; this excludes the Finish function's line template.
    used = set()
    for command in commands:
        used.update(re.findall(r"\$(\w+)", command))
    grew = True
    while grew:
        grew = False
        for name, value in values:
            if name in used:
                for inner in re.findall(r"\$(\w+)", value):
                    if inner not in used:
                        used.add(inner)
                        grew = True

    known = dict(values)
    def expand(value):
        for _ in range(5):
            wider = re.sub(
                r"\$(\w+)", lambda m: known.get(m.group(1), m.group(0)), value)
            if wider == value:
                return value
            value = wider
        return value

    all_ok = len(steps) == 9 and all(step[2] == 0 for step in steps)

    lines = [
        "Lane command summary (v5.2.0 bridge), "
        + datetime.now().astimezone().isoformat(timespec="seconds"),
        "",
        "Probe runs and analyses, FINISH lines as the lane script wrote them:",
    ]
    for _, step, code, seconds, name, line in steps:
        lines.append("  logs/%-26s %s" % (name, line))
    lines.append("")
    lines.append("%d steps, %s." % (
        len(steps), "all EXIT=0" if all_ok else "NOT all EXIT=0"))
    lines.append("")

    lines.append(
        "The lane script, %s, ran the commands below; the variables are "
        "resolved the way PowerShell resolves them:" % SCRIPT)
    for name, value in values:
        if name in used:
            lines.append("  $%s = '%s'" % (name, expand(value)))
    lines.append("")
    for command in commands:
        lines.append("  " + command)
    lines.append("")
    lines.append(
        "The six run logs hold the probe's console output and its FINISH line; "
        "each analysis step writes its console output to "
        "logs/analysis-<set>-v520.log and its FINISH line to "
        "logs/analysis-<set>.log.")
    lines.append("")
    lines.append(
        "The two probe binaries were built before this run; the build commands "
        "and logs are not in this directory.  compiled-tree-compare-v520.txt "
        "records what they were built from: their paths, sizes, SHA-256 digests "
        "and timestamps, the sources under external/boys compared file by file "
        "with the v5.2.0 tag's worktree, and the version each probe's own record "
        "states, 5.2.0 for all six arms.")

    text_out = "\n".join(lines) + "\n"
    if out is None:
        sys.stdout.write(text_out)
    else:
        out.write_text(text_out, encoding="utf-8", newline="\n")
    return 0 if all_ok else 1


if __name__ == "__main__":
    sys.exit(main())
