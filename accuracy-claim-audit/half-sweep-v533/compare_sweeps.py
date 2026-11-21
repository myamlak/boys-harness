"""Compare two half-sweep binaries cell by cell.

Both files were written by the same half_sweep.cpp, which enumerates arguments and orders in a
fixed order, so row i of one file is row i of the other. Reports per-format changed returned
bits, changed bounds, and the changed-cell identity at any cell named by a results file.
"""

import json
import struct
import sys

import numpy as np

ROW_DTYPE = np.dtype([
    ("format", "u1"),
    ("n", "u1"),
    ("argbits", "<u2"),
    ("gotbits", "<u2"),
    ("got", "<f8"),
    ("bound", "<f8"),
])
HEADER_SIZE = 20


def read(path):
    with open(path, "rb") as f:
        magic, version, nmax, rows = struct.unpack("<IIIQ", f.read(HEADER_SIZE))
        if magic != 0x50535748:
            raise SystemExit(f"{path}: bad magic 0x{magic:08x}")
        data = np.fromfile(f, dtype=ROW_DTYPE)
    if data.size != rows:
        raise SystemExit(f"{path}: header rows {rows} != file rows {data.size}")
    return {"version": version, "nmax": nmax, "rows": rows, "data": data}


def main():
    old_path, new_path, out_path = sys.argv[1], sys.argv[2], sys.argv[3]
    old = read(old_path)
    new = read(new_path)
    a, b = old["data"], new["data"]
    if a.size != b.size:
        raise SystemExit(f"row counts differ: {a.size} vs {b.size}")
    if not np.array_equal(a["format"], b["format"]) or not np.array_equal(a["n"], b["n"]) \
            or not np.array_equal(a["argbits"], b["argbits"]):
        raise SystemExit("cell order differs between the two sweeps")

    out = {
        "old": old_path,
        "new": new_path,
        "header": {"version_old": old["version"], "version_new": new["version"],
                   "nmax_old": old["nmax"], "nmax_new": new["nmax"],
                   "rows_old": int(old["rows"]), "rows_new": int(new["rows"])},
        "formats": {},
    }
    for fmt, name in ((0, "binary16"), (1, "bfloat16")):
        sel = a["format"] == fmt
        bits_changed = a["gotbits"][sel] != b["gotbits"][sel]
        bound_changed = a["bound"][sel] != b["bound"][sel]
        assert not np.any(a["argbits"] != b["argbits"])
        idx = np.flatnonzero(sel)[bits_changed]
        cells = [{"n": int(a["n"][i]), "arg_bits_hex": f"0x{int(a['argbits'][i]):04X}",
                  "old_got_bits_hex": f"0x{int(a['gotbits'][i]):04X}",
                  "new_got_bits_hex": f"0x{int(b['gotbits'][i]):04X}",
                  "old_got": float(a["got"][i]), "new_got": float(b["got"][i]),
                  "bound": float(b["bound"][i])}
                 for i in idx]
        out["formats"][name] = {
            "cells": int(sel.sum()),
            "changed_got_bits": int(bits_changed.sum()),
            "changed_bound": int(bound_changed.sum()),
            "changed_cells": cells,
        }

    # Mark whether each results-file worst/top cell changed.
    if len(sys.argv) > 4:
        out["results_cells"] = {}
        for res_path, label in ((sys.argv[4], "v511"), (sys.argv[5], "v533")):
            with open(res_path, "r", encoding="utf-8") as f:
                res = json.load(f)
            marks = []
            for name in ("binary16", "bfloat16"):
                for key in ("worst",):
                    w = res["formats"][name][key]
                    m = ((a["format"] == (0 if name == "binary16" else 1))
                         & (a["n"] == w["n"]) & (a["argbits"] == int(w["arg_bits_hex"], 16)))
                    i = int(np.flatnonzero(m)[0])
                    marks.append({
                        "results": label, "format": name, "slot": key,
                        "n": int(w["n"]), "arg_bits_hex": w["arg_bits_hex"],
                        "old_got_bits_hex": f"0x{int(a['gotbits'][i]):04X}",
                        "new_got_bits_hex": f"0x{int(b['gotbits'][i]):04X}",
                        "changed": bool(a["gotbits"][i] != b["gotbits"][i]),
                    })
                for k, t in enumerate(res["formats"][name].get("top_cells", [])):
                    m = ((a["format"] == (0 if name == "binary16" else 1))
                         & (a["n"] == t["n"]) & (a["argbits"] == int(t["arg_bits_hex"], 16)))
                    i = int(np.flatnonzero(m)[0])
                    marks.append({
                        "results": label, "format": name, "slot": f"top{k}",
                        "n": int(t["n"]), "arg_bits_hex": t["arg_bits_hex"],
                        "old_got_bits_hex": f"0x{int(a['gotbits'][i]):04X}",
                        "new_got_bits_hex": f"0x{int(b['gotbits'][i]):04X}",
                        "changed": bool(a["gotbits"][i] != b["gotbits"][i]),
                    })
            out["results_cells"][label] = marks

    with open(out_path, "w", encoding="utf-8") as f:
        json.dump(out, f, indent=2)
        f.write("\n")
    for name in ("binary16", "bfloat16"):
        r = out["formats"][name]
        print(f"{name}: {r['cells']} cells, {r['changed_got_bits']} changed returned bits, "
              f"{r['changed_bound']} changed bounds")


if __name__ == "__main__":
    main()
