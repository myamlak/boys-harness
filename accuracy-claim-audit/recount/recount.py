#!/usr/bin/env python3
"""Independent recount of the audit-pilot counts in Paper B (claims-as-rows).

The script reads only: it opens the pinned source trees and the audit's aggregate
record, and writes counts.json.  It never executes anything from the trees.

Definitions
-----------

B. Boost.Math 1.92.0 test-suite bounds.
   Scope: <boost>/libs/math/test/test_*.cpp.
   One bound = one textual ``add_expected_result(...)`` call with at least eight
   top-level comma-separated arguments (compiler, stdlib, platform, type, group,
   function, max eps, mean eps, and usually a description).  Comments are removed
   before parsing.  Preprocessor conditionals inside the argument list are
   resolved with BOOST_MATH_ENABLE_SYCL undefined and every other macro treated
   as undefined, which is the configuration the audit reports.  A robustness
   variant that keeps both branches and only drops the directive lines is also
   counted.  The rule counts textual calls, not the rows a test binary would run.

V. crlibm 1.0beta5 vector assertions.
   Scope: <crlibm>/tests/*.testdata.
   One vector = one non-comment data row whose first token is a rounding mode
   N, P, M or Z and which carries four hexadecimal words (x_hi x_lo y_hi y_lo)
   or six words (pow: x_hi x_lo y_hi y_lo z_hi z_lo).  A line whose first token
   is not a mode is a title or a header and is skipped; a mode line with any
   other word count or a non-hex word is reported as unparsed.  One file = one
   function.

G. crlibm Gappa obligations.
   Scope: <crlibm>/gappa/**/*.gappa.
   Primary count: one obligation = one line containing ``->``.  This is the
   rule that reproduces the paper's 791.  Alternatives are reported: arrow
   occurrences, non-comment arrow lines (``#`` starts a Gappa comment), and
   non-comment arrow lines that terminate the rule on the same line (end in
   ``;``).  Also counted: executable ``$GAPPA`` invocations in gappa/*.sh.
   A Gappa obligation in the strict sense is a rewrite rule or hint; 791 is the
   arrow-line count, not the number of goals Gappa would report.

C. crlibm numeric bounds in source comments.
   Scope: top-level <crlibm>/*.c and <crlibm>/*.h, i.e. the shipped function
   sources and their private headers.  tests/, gappa/, maple/, docs/ and
   scs_lib/ are excluded.
   One bound = one physical comment line containing both (a) a power-of-two
   numeric literal (2^n, 2^(n), 2**(-n), 2**n; the exponent may be decimal)
   and (b) one of the words error, accurate, accuracy, precision, bound
   (case-insensitive).  This is a line count: a repeated comment counts once
   per occurrence.
   The rule reproduces the paper's 41.  Two narrower readings are reported
   mechanically because the paper does not state the rule:
     * accuracy_word_only: (b) restricted to error/accurate/accuracy/precision
       (38 lines);
     * value_bound: accuracy_word_only minus the seven lines of the form
       "We must prove extra accuracy in the interval 0..sin(2^-18)", whose
       numeric literal is a domain endpoint rather than a bound value (31).
   The difference between 41 and 38 is three pow.c lines that bound a
   range-reduction exponent ("bounded by 2^11"), not an accuracy.

Registry and failure split
--------------------------
The 12 correct-rounding failures, their 4 sinh / 8 atanpi split, the 44
undecided reference points, the 361,286 queries and the 1,445,144 rounding-mode
results are read from the deposited aggregate counts.json.  The per-query file
m2_results.json that m2_claims.py writes is not in the deposited record, so the
split is checked against the aggregate only and is NOT recounted from per-query
rows.  counts.json records this distinction.

Usage
-----
  python -I recount.py \\
      --boost   C:/Users/myaml/dev/_boost-smoke/boost-1_92-src/boost_1_92_0 \\
      --boost-headers C:/Users/myaml/dev/_witness/third_party/boost_1_92_0 \\
      --boost-archive C:/Users/myaml/dev/_witness/downloads/boost_1_92_0.tar.gz \\
      --crlibm  C:/Users/myaml/dev/_witness/third_party/crlibm-1.0beta5 \\
      --crlibm-archive C:/Users/myaml/dev/_witness/downloads/crlibm-eb306379.tar.gz \\
      --registry <...>/claim-audit-crlibm/work/counts.json \\
      --out counts.json

--boost is the full Boost source tree (the witness tree is headers-only);
--boost-headers is optional and checks that the witness headers are identical.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import sys
import tarfile
from pathlib import Path

# --------------------------------------------------------------------------
# shared helpers
# --------------------------------------------------------------------------


def read_text(path: Path) -> str:
    return path.read_text(encoding="utf-8", errors="replace")


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def archive_provenance(path: Path) -> dict:
    """Hash an archive and read its top directory name without extracting it."""
    info = {
        "path": str(path),
        "bytes": path.stat().st_size,
        "sha256": sha256_file(path),
        "top_directory": None,
    }
    try:
        with tarfile.open(path, "r:gz") as tf:
            names = tf.getnames()
        if names:
            info["top_directory"] = names[0].split("/")[0]
    except (tarfile.TarError, OSError) as exc:
        info["error"] = f"{type(exc).__name__}: {exc}"
    return info


def compare_subtree(a: Path, b: Path, subdir: str) -> dict:
    """Compare every file under subdir in two trees by SHA-256."""
    base_a, base_b = a / subdir, b / subdir
    if not base_a.is_dir() or not base_b.is_dir():
        return {"comparable": False, "reason": f"missing {subdir} in one tree"}
    files_a = {
        p.relative_to(base_a).as_posix(): p for p in base_a.rglob("*") if p.is_file()
    }
    files_b = {
        p.relative_to(base_b).as_posix(): p for p in base_b.rglob("*") if p.is_file()
    }
    common = sorted(set(files_a) & set(files_b))
    identical = sum(
        1 for rel in common if sha256_file(files_a[rel]) == sha256_file(files_b[rel])
    )
    return {
        "comparable": True,
        "subdir": subdir,
        "files_in_a": len(files_a),
        "files_in_b": len(files_b),
        "files_compared": len(common),
        "identical": identical,
        "only_in_a": sorted(set(files_a) - set(files_b))[:20],
        "only_in_b": sorted(set(files_b) - set(files_a))[:20],
        "differing": [rel for rel in common if sha256_file(files_a[rel]) != sha256_file(files_b[rel])][:20],
    }


# --------------------------------------------------------------------------
# Boost.Math test-suite bounds (rule B)
# --------------------------------------------------------------------------

BOOST_VERSION_RE = re.compile(r"#define\s+BOOST_VERSION\s+(\d+)")
BOOST_CALL = re.compile(r"add_expected_result\s*\(")
# Macros known to the audit's reported configuration.  Everything else is
# treated as undefined, as the audit's parser does.
BOOST_MACROS = {"BOOST_MATH_ENABLE_SYCL": False}


def boost_version(root: Path) -> int | None:
    path = root / "boost" / "version.hpp"
    if not path.is_file():
        return None
    m = BOOST_VERSION_RE.search(read_text(path))
    return int(m.group(1)) if m else None


def strip_c_comments(text: str) -> str:
    """Remove C/C++ comments, preserving newlines so line numbers stay true."""
    out: list[str] = []
    i, n = 0, len(text)
    in_str = in_block = False
    while i < n:
        c = text[i]
        if in_str:
            out.append(c)
            if c == "\\" and i + 1 < n:
                out.append(text[i + 1])
                i += 2
                continue
            if c == '"':
                in_str = False
        elif in_block:
            if c == "*" and i + 1 < n and text[i + 1] == "/":
                in_block = False
                i += 2
                continue
            if c == "\n":
                out.append("\n")
        elif c == '"':
            in_str = True
            out.append(c)
        elif c == "/" and i + 1 < n and text[i + 1] == "/":
            j = text.find("\n", i)
            if j < 0:
                break
            i = j
            continue
        elif c == "/" and i + 1 < n and text[i + 1] == "*":
            in_block = True
            i += 2
            continue
        else:
            out.append(c)
        i += 1
    return "".join(out)


def split_call_args(text: str, start: int, resolve_conditionals: bool) -> tuple[list[str], int, list[str]]:
    """Return (top-level arguments, index after ')', preprocessor lines met).

    ``text[start:]`` begins just after the call's opening parenthesis.
    With resolve_conditionals=True the branch not taken is dropped; with False
    directive lines are treated as whitespace and both branches are kept.
    """
    args: list[str] = []
    cur: list[str] = []
    depth = 1
    i, n = start, len(text)
    in_str = False
    drop = False
    directives: list[str] = []
    while i < n:
        c = text[i]
        if not in_str and c == "#":
            j = text.find("\n", i)
            if j < 0:
                j = n
            line = text[i:j].strip()
            directives.append(line)
            if resolve_conditionals:
                body = line.lstrip("#").strip()
                parts = body.split()
                name = parts[0] if parts else ""
                if name in ("ifdef", "ifndef", "if"):
                    if name == "ifdef":
                        value = BOOST_MACROS.get(parts[1], False) if len(parts) > 1 else False
                    elif name == "ifndef":
                        value = not BOOST_MACROS.get(parts[1], False) if len(parts) > 1 else True
                    else:
                        value = True
                    drop = not value
                elif name in ("else", "elif"):
                    drop = not drop
                elif name == "endif":
                    drop = False
            i = j + 1 if j < n else j
            continue
        if drop:
            i += 1
            continue
        if in_str:
            cur.append(c)
            if c == "\\" and i + 1 < n:
                cur.append(text[i + 1])
                i += 2
                continue
            if c == '"':
                in_str = False
        elif c == '"':
            in_str = True
            cur.append(c)
        elif c == "(":
            depth += 1
            cur.append(c)
        elif c == ")":
            depth -= 1
            if depth == 0:
                args.append("".join(cur))
                return args, i + 1, directives
            cur.append(c)
        elif c == "," and depth == 1:
            args.append("".join(cur))
            cur = []
        else:
            cur.append(c)
        i += 1
    return args, i, directives


def count_boost(root: Path) -> dict:
    testdir = root / "libs" / "math" / "test"
    if not testdir.is_dir():
        raise SystemExit(
            f"{root}: no libs/math/test; pass --boost as a full Boost source tree "
            "(the witness tree under _witness/third_party is headers-only)"
        )
    files = sorted(testdir.glob("test_*.cpp"))
    records: list[list] = []
    directives: set[str] = set()
    plain = 0
    for path in files:
        text = strip_c_comments(read_text(path))
        for m in BOOST_CALL.finditer(text):
            args, _, dirs = split_call_args(text, m.end(), True)
            directives.update(dirs)
            if len(args) >= 8:
                line = text[: m.start()].count("\n") + 1
                records.append([path.name, line])
            args_plain, _, _ = split_call_args(text, m.end(), False)
            if len(args_plain) >= 8:
                plain += 1
    per_file: dict[str, int] = {}
    for name, _ in records:
        per_file[name] = per_file.get(name, 0) + 1
    return {
        "version_macro": boost_version(root),
        "test_directory": str(testdir),
        "test_files": len(files),
        "bound_records": len(records),
        "bound_records_plain_preprocessor_variant": plain,
        "directives_met": sorted(directives),
        "per_file": dict(sorted(per_file.items())),
        "records": records,
    }


# --------------------------------------------------------------------------
# crlibm vectors (rule V)
# --------------------------------------------------------------------------

CRLIBM_MODES = {"N", "P", "M", "Z"}
CRLIBM_WORD = re.compile(r"^(?:0x)?([0-9a-fA-F]{1,8})$")


def count_vectors(root: Path) -> dict:
    testdir = root / "tests"
    files = sorted(testdir.glob("*.testdata"))
    per_file: dict[str, dict] = {}
    unparsed: list[list] = []
    total = 0
    for path in files:
        n = 0
        by_mode: dict[str, int] = {}
        for lineno, raw in enumerate(read_text(path).splitlines(), 1):
            line = raw.split("#", 1)[0].strip()
            if not line:
                continue
            tokens = line.split()
            if tokens[0] not in CRLIBM_MODES:
                continue  # function title or header text, not a data row
            words = tokens[1:]
            if len(words) not in (4, 6) or any(not CRLIBM_WORD.match(w) for w in words):
                unparsed.append([path.name, lineno, line[:160]])
                continue
            n += 1
            by_mode[tokens[0]] = by_mode.get(tokens[0], 0) + 1
        per_file[path.name] = {"vectors": n, "by_mode": dict(sorted(by_mode.items()))}
        total += n
    return {
        "test_directory": str(testdir),
        "files": len(files),
        "total": total,
        "per_file": per_file,
        "unparsed": unparsed,
    }


# --------------------------------------------------------------------------
# crlibm Gappa arrows (rule G)
# --------------------------------------------------------------------------

def count_gappa(root: Path) -> dict:
    gdir = root / "gappa"
    files = sorted(gdir.rglob("*.gappa"))
    per_file: dict[str, dict] = {}
    arrow_lines = arrow_occurrences = noncomment = terminated = 0
    for path in files:
        lines = read_text(path).splitlines()
        a = occ = nc = term = 0
        for raw in lines:
            if "->" not in raw:
                continue
            a += 1
            occ += raw.count("->")
            if raw.lstrip().startswith("#"):
                continue
            nc += 1
            if raw.rstrip().endswith(";"):
                term += 1
        per_file[path.relative_to(gdir).as_posix()] = {
            "arrow_lines": a,
            "arrow_occurrences": occ,
            "noncomment_arrow_lines": nc,
            "noncomment_terminated": term,
        }
        arrow_lines += a
        arrow_occurrences += occ
        noncomment += nc
        terminated += term

    invocations: list[list] = []
    for sh in sorted(gdir.rglob("*.sh")):
        for lineno, raw in enumerate(read_text(sh).splitlines(), 1):
            line = raw.split("#", 1)[0]
            if not line.strip():
                continue
            if "$GAPPA" in line and not line.lstrip().startswith("GAPPA="):
                invocations.append([sh.relative_to(gdir).as_posix(), lineno, line.strip()])
    return {
        "gappa_directory": str(gdir),
        "files": len(files),
        "arrow_lines": arrow_lines,
        "arrow_occurrences": arrow_occurrences,
        "noncomment_arrow_lines": noncomment,
        "noncomment_terminated": terminated,
        "noncomment_continued": noncomment - terminated,
        "scripted_invocations": len(invocations),
        "invocations": invocations,
        "per_file": per_file,
    }


# --------------------------------------------------------------------------
# crlibm source-comment bounds (rule C)
# --------------------------------------------------------------------------

COMMENT_NUM2 = re.compile(
    r"2\s*\^\s*\(?\s*[+-]?\s*\d+(?:\.\d+)?\s*\)?|2\s*\*\*\s*[+-]?\s*\d+"
)
COMMENT_BOUND = re.compile(r"error|accurate|accuracy|precision|bound", re.I)
COMMENT_ACC = re.compile(r"error|accurate|accuracy|precision", re.I)
COMMENT_INTERVAL = re.compile(r"in\s+the\s+interval", re.I)


def comment_lines(path: Path):
    """Yield (line_no, text) for every physical C/C++ comment line."""
    src = read_text(path)
    i, n, line = 0, len(src), 1
    while i < n:
        c = src[i]
        if c == "\n":
            line += 1
            i += 1
        elif src.startswith("//", i):
            j = src.find("\n", i)
            if j < 0:
                j = n
            yield line, src[i + 2 : j]
            i = j
        elif src.startswith("/*", i):
            j = src.find("*/", i + 2)
            if j < 0:
                j = n
            body = src[i + 2 : j]
            start = line
            for k, text in enumerate(body.split("\n")):
                yield start + k, text
            line += body.count("\n")
            i = j + 2
        elif c == '"':
            j = i + 1
            while j < n:
                if src[j] == "\\":
                    j += 2
                    continue
                if src[j] == '"':
                    j += 1
                    break
                if src[j] == "\n":
                    line += 1
                j += 1
            i = j
        else:
            i += 1


def count_source_comments(root: Path) -> dict:
    records: list[list] = []
    for path in sorted(root.glob("*")):
        if path.suffix not in (".c", ".h"):
            continue
        for lineno, raw in comment_lines(path):
            text = raw.strip(" *\t")
            if not text:
                continue
            if COMMENT_NUM2.search(text) and COMMENT_BOUND.search(text):
                records.append([path.name, lineno, text])
    accuracy_only = [r for r in records if COMMENT_ACC.search(r[2])]
    interval = [r for r in accuracy_only if COMMENT_INTERVAL.search(r[2])]
    magnitude_only = [r for r in records if not COMMENT_ACC.search(r[2])]
    return {
        "scope": "top-level *.c and *.h; tests/, gappa/, maple/, docs/, scs_lib/ excluded",
        "lines": len(records),
        "accuracy_word_only": len(accuracy_only),
        "interval_statements": len(interval),
        "value_bound": len(accuracy_only) - len(interval),
        "magnitude_only": len(magnitude_only),
        "records": records,
    }


# --------------------------------------------------------------------------
# audit registry and paper comparison
# --------------------------------------------------------------------------

def parse_int(pattern: str, text: str) -> int | None:
    m = re.search(pattern, text)
    return int(m.group(1).replace(",", "")) if m else None


def read_registry(path: Path) -> dict:
    rec = json.loads(read_text(path))
    workdir = path.parent
    mcn = rec["machine_checked_not_in_manual"]
    m1 = rec["measurements"]["m1_crlibm_vs_own_vectors"]
    m2 = rec["measurements"]["m2_crlibm_vs_independent_reference"]
    validated = rec["reference"]["validated_against"]

    oracle = next((s for s in validated if "exact Fraction" in s), "")
    own = next((s for s in validated if "own tests" in s), "")

    data_files = sorted(
        p.name
        for p in workdir.iterdir()
        if p.is_file() and p.suffix.lower() in (".json", ".tsv", ".csv", ".log")
    )
    return {
        "path": str(path),
        "library": rec["library"]["name"],
        "version": rec["library"]["version"],
        "machine_checked_not_in_manual": {
            "vector_assertions": mcn["test_vector_assertions"],
            "vector_functions": mcn["test_vector_functions"],
            "gappa_obligations": mcn["gappa_goal_obligations"],
            "gappa_files": mcn["gappa_files"],
            "source_comment_bounds": mcn["numeric_bounds_in_source_comments"],
        },
        "documented_claims_total": rec["documented_claims"]["total"],
        "documented_internal_numeric_bound_sentences": rec["documented_claims"][
            "internal_numeric_bound_sentences"
        ],
        "m1": {
            "vectors": m1["vectors"],
            "disagreements": m1["disagreements"],
            "disagreements_on_finite_values": m1["disagreements_on_finite_values"],
            "note": m1.get("note", ""),
        },
        "m2": {
            "queries": m2["queries"],
            "rounding_mode_results_checked": m2["rounding_mode_results_checked"],
            "failures": m2["failures"],
            "sinh_failures": m2["sinh_failures"],
            "atanpi_failures": m2["atanpi_failures"],
            "undecided_reference_points": m2["undecided_reference_points"],
        },
        "table_s9_prose": {
            "oracle_mismatches": parse_int(r"(\d+)\s+mismatches", oracle),
            "oracle_trials": parse_int(r"in\s+([\d,]+)\s+trials", oracle),
            "oracle_overflow_boundary": parse_int(
                r"([\d,]+)\s+overflow-boundary", oracle
            ),
            "own_vector_disagreements": parse_int(r"(\d+)\s+disagreements", own),
            "own_vector_decided_approx": parse_int(
                r"\(~([\d,]+)\)", own
            ),
        },
        "defects": {
            "sinh": rec["defects"]["sinh"],
            "atanpi": rec["defects"]["atanpi"],
        },
        "per_query_results_deposited": (workdir / "m2_results.json").exists(),
        "m2_script_writes_m2_results": "m2_results.json"
        in read_text(workdir / "m2_claims.py")
        if (workdir / "m2_claims.py").is_file()
        else None,
        "conclusion_files_in_workdir": data_files,
    }


def build_comparison(boost: dict, crlibm: dict, registry: dict) -> list[dict]:
    v = crlibm["vectors"]
    g = crlibm["gappa"]
    c = crlibm["source_comments"]
    reg = registry["machine_checked_not_in_manual"]
    m1 = registry["m1"]
    m2 = registry["m2"]
    s9 = registry["table_s9_prose"]

    def row(claim, paper, recount, rule, location, basis="engine", note=""):
        return {
            "claim": claim,
            "paper": paper,
            "recount": recount,
            "agrees": paper == recount,
            "rule": rule,
            "location": location,
            "basis": basis,
            "note": note,
        }

    rows = [
        row(
            "boost.suite_bounds",
            499,
            boost["bound_records"],
            "B",
            "main.tex sec:problem-contrast; supplement S6",
            note="one add_expected_result call with >=8 top-level arguments",
        ),
        row(
            "boost.test_files",
            522,
            boost["test_files"],
            "B",
            "supplement S6 ('499 bounds across 522 test files')",
            note="count of libs/math/test/test_*.cpp",
        ),
        row(
            "crlibm.vector_assertions",
            63617,
            v["total"],
            "V",
            "main.tex sec:problem-contrast; supplement S6, Table S9",
        ),
        row(
            "crlibm.vector_functions",
            21,
            v["files"],
            "V",
            "supplement S6 ('over 21 functions')",
            note="one *.testdata file per function",
        ),
        row(
            "crlibm.gappa_obligations",
            791,
            g["arrow_lines"],
            "G",
            "main.tex sec:problem-contrast; supplement S6, S8",
            note=(
                "arrow-line count; alternatives: "
                f"{g['noncomment_arrow_lines']} non-comment arrow lines, "
                f"{g['noncomment_terminated']} arrow lines ending in ';', "
                f"{g['files']} files, {g['scripted_invocations']} scripted gappa invocations"
            ),
        ),
        row(
            "crlibm.gappa_files",
            46,
            g["files"],
            "G",
            "supplement S6 ('across 46 files')",
        ),
        row(
            "crlibm.source_comment_bounds",
            41,
            c["lines"],
            "C",
            "main.tex sec:problem-contrast; supplement S6",
            note=(
                "the broad line rule reproduces 41; the paper does not state it. "
                f"accuracy-word-only reading: {c['accuracy_word_only']}; "
                f"value-bound reading (accuracy word minus {c['interval_statements']} "
                f"interval statements): {c['value_bound']}; "
                f"{c['magnitude_only']} magnitude-only lines ('bounded by 2^11')"
            ),
        ),
        row(
            "crlibm.correct_rounding_failures",
            12,
            m2["failures"],
            "registry",
            "supplement Table S9; main.tex sec:problem-held",
            basis="registry aggregate only",
            note="per-query results not deposited",
        ),
        row(
            "crlibm.failures.sinh",
            4,
            m2["sinh_failures"],
            "registry",
            "supplement Table S9; main.tex sec:problem-held",
            basis="registry aggregate only",
            note="per-query results not deposited",
        ),
        row(
            "crlibm.failures.atanpi",
            8,
            m2["atanpi_failures"],
            "registry",
            "supplement Table S9; main.tex sec:problem-held",
            basis="registry aggregate only",
            note="registry describes these as outside the documented domain",
        ),
        row(
            "crlibm.undecided_reference_points",
            44,
            m2["undecided_reference_points"],
            "registry",
            "supplement Table S9",
            basis="registry aggregate only",
            note="m2 field; not the unrelated 44 internal manual sentences",
        ),
        row(
            "crlibm.queries",
            361286,
            m2["queries"],
            "registry",
            "supplement Table S9",
            basis="registry aggregate only",
            note=f"4 rounding modes per query: {m2['queries']} * 4 = {m2['rounding_mode_results_checked']}",
        ),
        row(
            "crlibm.rounding_mode_results",
            1445144,
            m2["rounding_mode_results_checked"],
            "registry",
            "supplement Table S9",
            basis="registry aggregate only",
        ),
        row(
            "table_s9.oracle_mismatches",
            0,
            s9["oracle_mismatches"],
            "registry",
            "supplement Table S9",
            basis="registry prose",
        ),
        row(
            "table_s9.oracle_trials",
            5999,
            s9["oracle_trials"],
            "registry",
            "supplement Table S9",
            basis="registry prose",
        ),
        row(
            "table_s9.oracle_overflow_boundary",
            72,
            s9["oracle_overflow_boundary"],
            "registry",
            "supplement Table S9",
            basis="registry prose",
        ),
        row(
            "table_s9.own_vectors_decided_approx",
            48000,
            s9["own_vector_decided_approx"],
            "registry",
            "supplement Table S9 ('about 48,000')",
            basis="registry prose; stated as approximate",
            note="the paper says about 48,000; the registry says (~48,000)",
        ),
        row(
            "table_s9.own_vector_disagreements",
            0,
            s9["own_vector_disagreements"],
            "registry",
            "supplement Table S9",
            basis="registry prose",
        ),
        row(
            "table_s9.vector_count",
            63617,
            v["total"],
            "V",
            "supplement Table S9",
            note="the registry m1 vectors field is also 63,617",
        ),
        row(
            "table_s9.vector_disagreements",
            9,
            m1["disagreements"],
            "registry",
            "supplement Table S9",
            basis="registry aggregate only",
            note=f"{m1['disagreements_on_finite_values']} on finite values; all 9 are NaN bit patterns",
        ),
    ]
    return rows


# --------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------

def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--boost", required=True, type=Path, help="full Boost source tree")
    ap.add_argument("--boost-headers", type=Path, help="pinned witness headers tree")
    ap.add_argument("--boost-archive", type=Path)
    ap.add_argument("--boost-sha256", default="c4a3b310ddd2472416e091067166b0713be97c63f38c212c484ada022fd296ce")
    ap.add_argument("--crlibm", required=True, type=Path, help="crlibm source tree")
    ap.add_argument("--crlibm-archive", type=Path)
    ap.add_argument("--crlibm-commit", default="eb3063791aa75bc9705b49283bf14250465220a7")
    ap.add_argument("--registry", required=True, type=Path, help="audit counts.json")
    ap.add_argument("--out", type=Path, default=Path("counts.json"))
    args = ap.parse_args(argv)

    # provenance
    provenance: dict = {}
    bv = boost_version(args.boost)
    provenance["boost"] = {
        "source_root": str(args.boost),
        "version_macro": bv,
        "version_ok": bv == 109200,
        "has_libs_math_test": (args.boost / "libs" / "math" / "test").is_dir(),
    }
    if args.boost_archive:
        prov = archive_provenance(args.boost_archive)
        prov["expected_sha256"] = args.boost_sha256
        prov["sha256_matches"] = prov["sha256"] == args.boost_sha256
        provenance["boost"]["archive"] = prov
    if args.boost_headers:
        provenance["boost"]["witness_headers"] = {
            "root": str(args.boost_headers),
            "version_macro": boost_version(args.boost_headers),
            "boost_math": compare_subtree(args.boost, args.boost_headers, "boost/math"),
            "has_libs_math_test": (args.boost_headers / "libs" / "math" / "test").is_dir(),
        }

    version_file = args.crlibm / "VERSION"
    cr_version = read_text(version_file).strip() if version_file.is_file() else None
    provenance["crlibm"] = {
        "source_root": str(args.crlibm),
        "version_file": cr_version,
        "version_ok": cr_version == "1.0beta5",
        "expected_commit": args.crlibm_commit,
    }
    if args.crlibm_archive:
        prov = archive_provenance(args.crlibm_archive)
        top = prov.get("top_directory") or ""
        prov["commit_from_top_directory"] = top.split("-", 1)[1] if top.startswith("crlibm-") else None
        prov["commit_matches"] = prov["commit_from_top_directory"] == args.crlibm_commit
        provenance["crlibm"]["archive"] = prov

    # counts
    boost = count_boost(args.boost)
    boost["version_ok"] = bv == 109200
    crlibm = {
        "vectors": count_vectors(args.crlibm),
        "gappa": count_gappa(args.crlibm),
        "source_comments": count_source_comments(args.crlibm),
    }
    registry = read_registry(args.registry)
    comparison = build_comparison(boost, crlibm, registry)

    recalculable = [r for r in comparison if r["basis"] == "engine"]
    aggregate_only = [r for r in comparison if r["basis"].startswith("registry")]
    summary = {
        "claims_checked": len(comparison),
        "engine_recounted": len(recalculable),
        "engine_agree": sum(1 for r in recalculable if r["agrees"]),
        "aggregate_only": len(aggregate_only),
        "aggregate_agree": sum(1 for r in aggregate_only if r["agrees"]),
        "disagreements": [r["claim"] for r in comparison if not r["agrees"]],
        "per_query_failure_split_recounted": False,
    }

    out = {
        "paper": "papers-b-v5/papers/claims-as-rows",
        "generated_by": "recount.py",
        "definitions": {
            "B": "add_expected_result(...) call with >=8 top-level arguments in libs/math/test/test_*.cpp",
            "V": "data row in tests/*.testdata with mode N/P/M/Z and 4 or 6 hex words",
            "G": "line containing '->' in gappa/**/*.gappa (alternatives reported)",
            "C": "comment line in top-level *.c/*.h with a power-of-two literal and one of error/accurate/accuracy/precision/bound",
        },
        "provenance": provenance,
        "boost": boost,
        "crlibm": crlibm,
        "registry": registry,
        "comparison": comparison,
        "summary": summary,
    }
    args.out.write_text(json.dumps(out, indent=2) + "\n", encoding="utf-8")

    print(f"wrote {args.out}")
    print(
        f"Boost {bv}: {boost['test_files']} test files, "
        f"{boost['bound_records']} bounds (plain variant {boost['bound_records_plain_preprocessor_variant']})"
    )
    print(
        f"crlibm {cr_version}: {crlibm['vectors']['total']} vectors over "
        f"{crlibm['vectors']['files']} files; Gappa {crlibm['gappa']['arrow_lines']} arrow lines "
        f"over {crlibm['gappa']['files']} files; source comments {crlibm['source_comments']['lines']}"
    )
    print(
        f"registry: {registry['m2']['failures']} failures "
        f"({registry['m2']['sinh_failures']} sinh + {registry['m2']['atanpi_failures']} atanpi), "
        f"{registry['m2']['undecided_reference_points']} undecided, "
        f"{registry['m2']['queries']} queries; per-query results deposited: "
        f"{registry['per_query_results_deposited']}"
    )
    print(
        f"comparison: {summary['engine_recounted']} engine-recounted "
        f"({summary['engine_agree']} agree), {summary['aggregate_only']} aggregate-only "
        f"({summary['aggregate_agree']} agree), disagreements: {summary['disagreements'] or 'none'}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
