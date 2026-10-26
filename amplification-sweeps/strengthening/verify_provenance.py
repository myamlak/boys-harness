"""Provenance checks for the strengthening analysis, from the deposit alone.

    python verify_provenance.py

The original provenance run, recorded byte for byte in verify_provenance.log,
also compared the paper tree's copy of analyze.py with the deposited copy.  A
clone of this deposit has no paper tree, so that half cannot be repeated here;
the committed log remains the record of the original run.  This deposit-only
variant checks, each printed with PASS or FAIL and reflected in the exit code:

  1. the three decompressed inputs under inputs/ are byte-identical to the
     deposit's xz streams (the xz streams themselves are never modified);
  2. the deposited analyze.py has the sha256 the original run printed;
  3. the two provenance.json records name the same boys commit, tag, compiler,
     compile flags (except the source and output paths) and the same 15 boys
     header-test-tool digests.

Its output therefore differs from the committed verify_provenance.log in two
ways: the three input labels say inputs/ rather than raw/, and the three
paper-tree body-identity PASS lines and the paper-tree hash line are absent.
No hash or provenance value changes.
"""

import hashlib
import json
import lzma
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
DEPOSIT = HERE.parent
ANALYZE_SHA256 = "4e7454492a53cef13b7ce2b98ff10bde0529ac406c105a570193499a69077f20"
failures = []


def check(label, ok, detail=""):
    print("%-4s %s%s" % ("PASS" if ok else "FAIL", label, ("  " + detail) if detail else ""))
    if not ok:
        failures.append(label)


def stream_hash(handle):
    digest = hashlib.sha256()
    while True:
        block = handle.read(1 << 22)
        if not block:
            return digest.hexdigest()
        digest.update(block)


def file_hash(path):
    with open(path, "rb") as handle:
        return stream_hash(handle)


def main():
    # ---- 1. inputs against the deposit's streams --------------------------
    pairs = (
        ("inputs/cells.csv", "product-amplification/results/cells.csv.xz"),
        ("inputs/reference.csv", "product-amplification/data/reference.csv.xz"),
        ("inputs/paths.csv", "amplification-form/results/paths.csv.xz"),
    )
    for local, packed in pairs:
        if not (HERE / local).is_file():
            check("%s == %s" % (local, packed), False, "missing; run prepare_inputs.py")
            continue
        local_hash = file_hash(HERE / local)
        with lzma.open(DEPOSIT / packed, "rb") as handle:
            deposit_hash = stream_hash(handle)
        check("%s == %s" % (local, packed), local_hash == deposit_hash,
              local_hash if local_hash == deposit_hash else
              "%s vs %s" % (local_hash, deposit_hash))

    # ---- 2. the deposited analyze.py --------------------------------------
    analyze_hash = file_hash(DEPOSIT / "amplification-form" / "analyze.py")
    check("analyze.py sha256", analyze_hash == ANALYZE_SHA256, analyze_hash)

    # ---- 3. the two provenance records ------------------------------------
    first = json.load(open(DEPOSIT / "product-amplification" / "results" / "provenance.json",
                           encoding="utf-8"))
    second = json.load(open(DEPOSIT / "amplification-form" / "results" / "provenance.json",
                            encoding="utf-8"))
    check("same boys commit", first["boys"]["commit"] == second["boys"]["commit"],
          first["boys"]["commit"])
    check("same boys tag", first["boys"]["tag_at_head"] == second["boys"]["tag_at_head"],
          first["boys"]["tag_at_head"])
    check("both boys trees clean",
          not first["boys"]["dirty"] and not second["boys"]["dirty"])
    check("same compiler", first["compiler"]["path"] == second["compiler"]["path"],
          first["compiler"]["path"])
    flags_first = first["flags"].split(" /I ")[0]
    flags_second = second["flags"].split(" /I ")[0]
    check("same compile flags", flags_first == flags_second, flags_first)
    keys_first = {k: v for k, v in first["digests"].items() if k.startswith("boys/")}
    keys_second = {k: v for k, v in second["digests"].items() if k.startswith("boys/")}
    check("same 15 boys digests", keys_first == keys_second,
          "%d entries" % len(keys_first))
    check("same machine",
          first["machine"]["node"] == second["machine"]["node"]
          and first["machine"]["processor"] == second["machine"]["processor"],
          first["machine"]["node"])

    print()
    print("all checks passed" if not failures else "FAILURES: %s" % failures)
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
