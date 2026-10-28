#!/usr/bin/env python3
"""Check the documented-claims enumeration against the gate's registration table.

    python accuracy-claim-audit/tools/check_registration.py --boys <boys-checkout>

The gate source is read at <boys-checkout>/tests/boys_accuracy_gate.cpp and
the enumeration defaults to records/documented-claims-v5.tsv beside this
directory; --tsv names another copy.

The checker parses the claims book's add() call sites in printed order (the
enumeration's G indices), the registration table kRegistrations[], and the
book's row identifiers, under the preprocessor configuration the published
numbering was made in: BOYS_GATE_FP16 defined, BOYS_GATE_EXTERNAL_BOOST not
defined. It prints the G-index mapping it used.

It checks:

(a) the set of claim ids with status claims-book in the enumeration equals
    the set of claim ids in the registration table;
(b) for every claims-book claim, every registered_by G index is a book row,
    the row the gate registers for the claim is cited by the claim, and every
    cited row is registered to the claim, unregistered, or shared with
    another claim that cites the same row. A registration to a row appended
    after the enumeration's highest cited index cannot be cited by the frozen
    enumeration; it is accepted only when each cited row is marked in
    kUnregisteredRows or cited by another claim registered to it, and is
    printed as a supersession notice, or counted as a mismatch under --strict;
(c) every row named by a registration is a row of the claims book.

Exit status 0 when the two agree, 2 on any mismatch or parse failure.

Only the Python standard library is used.
"""

import argparse
import csv
import re
import sys
from pathlib import Path

GATE_RELATIVE = Path("tests") / "boys_accuracy_gate.cpp"
DEFAULT_TSV = Path(__file__).resolve().parent.parent / "records" / "documented-claims-v5.tsv"

# The configuration the published claims-book numbering was made in.
DEFINES = {
    "BOYS_GATE_FP16": True,
    "BOYS_GATE_EXTERNAL_BOOST": False,
    "__x86_64__": True,
    "_M_X64": True,
}

CONDITIONAL = re.compile(r"^\s*#\s*(ifdef|ifndef|if|elif|else|endif)\b(.*)$")
DEFINED_ATOM = re.compile(r"defined\s*\(\s*([A-Za-z_][A-Za-z0-9_]*)\s*\)$")
G_TOKEN = re.compile(r"^G([0-9]+)$")
REGISTRATION_ENTRY = re.compile(r'\{\s*"([^"]*)"\s*,\s*"([^"]*)"\s*,')
BOOK_DECLARATION = re.compile(r"std::vector<\s*DocClaim\s*>\s+book\s*;")
REGISTRATION_STRUCT = re.compile(r"\bstruct\s+Registration\b")
REGISTRATION_ARRAY = re.compile(r"\bkRegistrations\s*\[\s*\]\s*=")
UNREGISTERED_ARRAY = re.compile(r"\bkUnregisteredRows\s*\[\s*\]\s*=")
UNREGISTERED_ENTRY = re.compile(r'\{\s*"([^"]+)"\s*,')


class ParseError(Exception):
    """A source or enumeration that cannot be read as the checker expects."""


class NoRegistrationTable(ParseError):
    """The gate checkout has no registration table."""


def _atom(text):
    text = text.strip()
    match = DEFINED_ATOM.match(text)
    if match:
        return bool(DEFINES.get(match.group(1), False))
    if text == "1":
        return True
    if text == "0":
        return False
    return bool(DEFINES.get(text, False))


def _condition(kind, expression):
    expression = expression.strip()
    if kind == "ifdef":
        name = expression.split()[0] if expression else ""
        return bool(DEFINES.get(name, False))
    if kind == "ifndef":
        name = expression.split()[0] if expression else ""
        return not bool(DEFINES.get(name, False))
    if kind == "if":
        if "||" in expression:
            return any(_atom(part) for part in expression.split("||"))
        if "&&" in expression:
            return all(_atom(part) for part in expression.split("&&"))
        return _atom(expression)
    raise ParseError("unknown conditional directive")


def active_lines(lines):
    """Yield (line number, line) for the lines that execute under DEFINES."""
    stack = []
    active = True
    for number, line in enumerate(lines, 1):
        match = CONDITIONAL.match(line)
        if match:
            kind, expression = match.group(1), match.group(2)
            if kind in ("ifdef", "ifndef", "if"):
                condition = _condition(kind, expression)
                stack.append([active, condition, active and condition])
                active = active and condition
            elif kind == "elif":
                if not stack:
                    raise ParseError("a conditional has no opening directive")
                parent, taken, _ = stack[-1]
                condition = (not taken) and _condition("if", expression)
                stack[-1][1] = taken or condition
                stack[-1][2] = parent and condition
                active = parent and condition
            elif kind == "else":
                if not stack:
                    raise ParseError("a conditional has no opening directive")
                parent, taken, _ = stack[-1]
                condition = not taken
                stack[-1][1] = True
                stack[-1][2] = parent and condition
                active = parent and condition
            else:
                if not stack:
                    raise ParseError("an #endif has no opening directive")
                parent, _, _ = stack.pop()
                active = parent
            continue
        if active:
            yield number, line
    if stack:
        raise ParseError("a conditional is not closed")


def _string_literals(text):
    """Return (start, end, value) for every C++ string literal in text."""
    literals = []
    index = 0
    while index < len(text):
        if text[index] != '"':
            index += 1
            continue
        cursor = index + 1
        chars = []
        while cursor < len(text):
            char = text[cursor]
            if char == "\\":
                following = text[cursor + 1] if cursor + 1 < len(text) else ""
                chars.append(
                    {"n": "\n", "t": "\t", "r": "\r", "0": "\0"}.get(following, following)
                )
                cursor += 2
                continue
            if char == '"':
                break
            chars.append(char)
            cursor += 1
        literals.append((index, cursor + 1, "".join(chars)))
        index = cursor + 1
    return literals


def first_string_argument(text):
    """The concatenated string literal in the first argument, or None."""
    open_paren = text.find("(")
    if open_paren < 0:
        return None
    depth = 0
    index = open_paren + 1
    while index < len(text):
        char = text[index]
        if char == '"':
            cursor = index + 1
            while cursor < len(text):
                if text[cursor] == "\\":
                    cursor += 2
                    continue
                if text[cursor] == '"':
                    break
                cursor += 1
            index = cursor + 1
            continue
        if char == "(":
            depth += 1
        elif char == ")":
            if depth == 0:
                break
            depth -= 1
        elif char == "," and depth == 0:
            break
        index += 1
    pieces = [value for _, _, value in _string_literals(text[open_paren + 1 : index])]
    return "".join(pieces) if pieces else None


def read_gate(path):
    if not path.is_file():
        raise ParseError(f"{path}: the gate source does not exist")
    try:
        return path.read_text(encoding="utf-8")
    except OSError as exc:
        raise ParseError(f"{path}: cannot be read: {exc}") from exc
    except UnicodeDecodeError as exc:
        raise ParseError(f"{path}: is not UTF-8 text: {exc}") from exc


def parse_registrations(text, path):
    if not REGISTRATION_STRUCT.search(text) and not REGISTRATION_ARRAY.search(text):
        raise NoRegistrationTable(
            f"{path}: no registration table (struct Registration / kRegistrations[]); "
            "this gate carries none"
        )
    if not REGISTRATION_STRUCT.search(text) or not REGISTRATION_ARRAY.search(text):
        raise ParseError(f"{path}: the registration table is incomplete")
    lines = text.splitlines()
    start = None
    for number, line in enumerate(lines, 1):
        if REGISTRATION_ARRAY.search(line):
            start = number
            break
    if start is None:
        raise ParseError(f"{path}: kRegistrations[] is not in the source")
    end = None
    for number in range(start + 1, len(lines) + 1):
        if lines[number - 1].strip() == "};":
            end = number
            break
    if end is None:
        raise ParseError(f"{path}: the registration table has no closing }};")
    region = "\n".join(lines[start - 1 : end])
    entries = []
    for match in REGISTRATION_ENTRY.finditer(region):
        claim, row = match.group(1), match.group(2)
        if claim and row:
            entries.append((claim, row))
    if not entries:
        raise ParseError(f"{path}: the registration table has no entries")
    return entries


def parse_unregistered(text, path):
    """The rows the gate explicitly marks as carrying no registered claim."""
    lines = text.splitlines()
    start = None
    for number, line in enumerate(lines, 1):
        if UNREGISTERED_ARRAY.search(line):
            start = number
            break
    if start is None:
        return set()
    end = None
    for number in range(start + 1, len(lines) + 1):
        if lines[number - 1].strip() == "};":
            end = number
            break
    if end is None:
        raise ParseError(f"{path}: kUnregisteredRows[] has no closing }};")
    region = "\n".join(lines[start - 1 : end])
    return {match.group(1) for match in UNREGISTERED_ENTRY.finditer(region)}


def parse_book(text, path):
    lines = text.splitlines()
    start = None
    end = None
    for number, line in enumerate(lines, 1):
        if start is None and BOOK_DECLARATION.search(line):
            start = number + 1
        if start is not None and REGISTRATION_STRUCT.search(line):
            end = number
            break
    if start is None:
        raise ParseError(
            f"{path}: the claims book (std::vector<DocClaim> book;) is not in the source"
        )
    if end is None:
        raise ParseError(f"{path}: the claims book has no registration table after it")
    rows = []
    seen = {}
    for number, line in active_lines(lines):
        if number < start or number >= end:
            continue
        stripped = line.strip()
        if not (stripped.startswith("add(") or stripped.startswith("add (")):
            continue
        row = first_string_argument(stripped)
        if row is None:
            continue
        if row in seen:
            raise ParseError(
                f'{path}: the claims book adds the row "{row}" twice '
                f"(lines {seen[row]} and {number})"
            )
        seen[row] = number
        rows.append((row, number))
    if not rows:
        raise ParseError(f"{path}: no claims-book add() call sites found")
    return rows


def read_tsv(path):
    try:
        with path.open(encoding="utf-8-sig", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if reader.fieldnames is None:
                raise ParseError(f"{path}: the enumeration is empty")
            required = ("id", "status", "registered_by")
            missing = [column for column in required if column not in reader.fieldnames]
            if missing:
                raise ParseError(
                    f"{path}: the enumeration has no column(s): {', '.join(missing)}"
                )
            rows = [row for row in reader if any(value for value in row.values() if value)]
    except OSError as exc:
        raise ParseError(f"{path}: cannot be read: {exc}") from exc
    except csv.Error as exc:
        raise ParseError(f"{path}: is not a readable TSV: {exc}") from exc
    if not rows:
        raise ParseError(f"{path}: the enumeration has no data rows")
    status_by_id = {}
    claims = []
    for row in rows:
        claim = (row.get("id") or "").strip()
        status = (row.get("status") or "").strip()
        if not claim:
            raise ParseError(f"{path}: a data row has no id")
        if claim in status_by_id:
            raise ParseError(f'{path}: the claim id "{claim}" appears twice')
        status_by_id[claim] = status
        if status == "claims-book":
            claims.append(row)
    return claims, status_by_id


def parse_g_labels(text, claim, path):
    labels = []
    raw = (text or "").strip()
    if not raw:
        return labels
    for token in re.split(r"[;,]", raw):
        token = token.strip()
        if not token:
            continue
        match = G_TOKEN.match(token)
        if not match:
            raise ParseError(
                f'{path}: claim {claim}: registered_by entry "{token}" is not a G index'
            )
        index = int(match.group(1))
        if index < 1:
            raise ParseError(
                f"{path}: claim {claim}: registered_by entry {token} is not a G index"
            )
        labels.append((index, token))
    return labels


def cite_text(labels, g_to_row):
    parts = []
    for index, token in labels:
        row = g_to_row.get(index)
        parts.append(f'{token}="{row}"' if row is not None else token)
    return ", ".join(parts)


def registration_text(rows, row_to_g):
    parts = []
    for row in rows:
        index = row_to_g.get(row)
        parts.append(
            f'"{row}" (G{index})' if index is not None else f'"{row}" (not in the book)'
        )
    return ", ".join(parts)


def main():
    parser = argparse.ArgumentParser(
        description="Check the documented-claims enumeration against the gate's registrations."
    )
    parser.add_argument(
        "--boys",
        required=True,
        type=Path,
        help="the boys checkout holding tests/boys_accuracy_gate.cpp",
    )
    parser.add_argument(
        "--tsv",
        type=Path,
        default=DEFAULT_TSV,
        help="the documented-claims enumeration (default: records/documented-claims-v5.tsv)",
    )
    parser.add_argument(
        "--strict",
        action="store_true",
        help="count a registration to a row appended after the enumeration's reach as a mismatch",
    )
    args = parser.parse_args()

    gate = args.boys / GATE_RELATIVE
    try:
        text = read_gate(gate)
        registrations = parse_registrations(text, gate)
        released = parse_unregistered(text, gate)
        book = parse_book(text, gate)
        claims, status_by_id = read_tsv(args.tsv)

        g_to_row = {index: row for index, (row, _) in enumerate(book, 1)}
        row_to_g = {row: index for index, row in g_to_row.items()}

        registered = {}
        owners = {}
        for claim, row in registrations:
            registered.setdefault(claim, [])
            if row not in registered[claim]:
                registered[claim].append(row)
            owners.setdefault(row, [])
            if claim not in owners[row]:
                owners[row].append(claim)

        claim_ids = [row["id"].strip() for row in claims]
        claim_ids_set = set(claim_ids)
        tsv_registered_by = {
            row["id"].strip(): (row.get("registered_by") or "") for row in claims
        }

        mismatches = []
        notices = []

        # (a) the two claim sets
        for row in claims:
            claim = row["id"].strip()
            if claim not in registered:
                mismatches.append(
                    f"claim {claim}: TSV entry status=claims-book, "
                    f'registered_by="{tsv_registered_by[claim]}"; '
                    "gate entry: no registration"
                )
        for claim in registered:
            if claim not in claim_ids_set:
                status = status_by_id.get(claim, "absent")
                mismatches.append(
                    f"claim {claim}: TSV entry status={status}; "
                    f"gate entry: {registration_text(registered[claim], row_to_g)}"
                )

        # (b) the citations
        cited = {}
        highest = 0
        for row in claims:
            claim = row["id"].strip()
            labels = parse_g_labels(row.get("registered_by"), claim, args.tsv)
            cited[claim] = labels
            for index, _ in labels:
                if 1 <= index <= len(book):
                    highest = max(highest, index)
        reach = highest if highest else len(book)

        for row in claims:
            claim = row["id"].strip()
            labels = cited[claim]
            if claim not in registered:
                continue
            if not labels:
                mismatches.append(
                    f'claim {claim}: TSV entry registered_by="" (no G index); '
                    f"gate entry: {registration_text(registered[claim], row_to_g)}"
                )
                continue
            cited_rows = {}
            for index, token in labels:
                if not (1 <= index <= len(book)):
                    mismatches.append(
                        f"claim {claim}: TSV entry registered_by="
                        f'"{tsv_registered_by[claim]}"; '
                        f"gate entry: {token} is beyond the book's {len(book)} rows"
                    )
                    continue
                cited_rows[index] = g_to_row[index]
            for registered_row in registered[claim]:
                if registered_row not in row_to_g:
                    continue  # (c) reports it
                index = row_to_g[registered_row]
                if index > reach:
                    stale = []
                    for cited_index, cited_row in sorted(cited_rows.items()):
                        if cited_row in released:
                            continue
                        row_owners = owners.get(cited_row, [])
                        if any(
                            cited_index in {i for i, _ in cited.get(owner, [])}
                            for owner in row_owners
                            if owner != claim
                        ):
                            continue
                        stale.append(f'G{cited_index}="{cited_row}"')
                    if stale:
                        mismatches.append(
                            f"claim {claim}: TSV entry cites {', '.join(stale)}; "
                            f'gate entry: registered to "{registered_row}" (G{index}), '
                            f"appended beyond G{reach}, and the cited row is neither marked "
                            "in kUnregisteredRows nor cited by another claim registered to it"
                        )
                        continue
                    note = (
                        f"claim {claim}: TSV entry cites {cite_text(labels, g_to_row)}; "
                        f'gate entry: registered to "{registered_row}" (G{index}), appended '
                        f"beyond the enumeration's highest cited index G{reach}"
                    )
                    if args.strict:
                        mismatches.append(note)
                    else:
                        notices.append(note)
                    continue
                if index not in cited_rows:
                    mismatches.append(
                        f"claim {claim}: TSV entry cites {cite_text(labels, g_to_row)}; "
                        f'gate entry: registered to "{registered_row}" (G{index}), '
                        "which the claim does not cite"
                    )
            for index, cited_row in sorted(cited_rows.items()):
                row_owners = owners.get(cited_row, [])
                if claim in row_owners or not row_owners:
                    continue
                shared = any(
                    index in {i for i, _ in cited.get(owner, [])} for owner in row_owners
                )
                if shared:
                    continue
                mismatches.append(
                    f'claim {claim}: TSV entry cites G{index}="{cited_row}"; '
                    f"gate entry: that row is registered to {', '.join(row_owners)}, "
                    f"which does not cite G{index}"
                )

        # (c) every registered row exists in the book
        for claim, row in registrations:
            if row not in row_to_g:
                mismatches.append(
                    f"claim {claim}: TSV entry registered_by="
                    f'"{tsv_registered_by.get(claim, "")}"; '
                    f'gate entry: registration row "{row}" is not a row of the claims book'
                )
    except ParseError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2

    reach_text = f"G{reach}" if highest else "no valid G index"
    print(f"gate source  : {gate}")
    print(f"enumeration  : {args.tsv}")
    print(f"book         : {len(book)} rows; the enumeration reaches {reach_text}")
    print(f"registrations: {len(registrations)} entries over {len(registered)} claims")
    print(f"claims-book  : {len(claims)} claims")
    print()
    print("the G mapping used (the claims book's add() order):")
    for index, (row, line) in enumerate(book, 1):
        print(f"  G{index:<3} {row:<42} line {line}")
    print()
    for notice in notices:
        print(f"notice:   {notice}")
    for mismatch in mismatches:
        print(f"mismatch: {mismatch}")

    if mismatches:
        print(f"\nregistration check: {len(mismatches)} mismatch(es)")
        return 2
    print(
        f"\nregistration check: OK ({len(claims)} claims-book claims, "
        f"{len(registrations)} registrations)"
    )
    if notices:
        print(f"  {len(notices)} superseded registration(s) noted")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
