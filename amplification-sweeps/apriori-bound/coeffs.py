"""Parse the v5.0.0 double-lane region-A tables straight from the headers.

The numbers are taken from the header text and converted through a Python
float, so the value used here is the exact binary double the compiler stores,
not the 18-digit decimal as written.
"""

import re

BOYS = r"C:/Users/myaml/dev/bladawiec/external/boys/include/boys"


def _read(name):
    with open(f"{BOYS}/{name}", "r", encoding="utf-8") as handle:
        return handle.read()


def _array_body(text, name):
    """The text between `name = std::to_array<...>({` and the closing `});`."""
    start = text.index(name)
    open_token = text.index("({", start)
    depth = 0
    for i in range(open_token + 1, len(text)):
        ch = text[i]
        if ch in "({":
            depth += 1
        elif ch in ")}":
            depth -= 1
            if depth == 0:
                return text[open_token + 2 : i]
    raise ValueError(f"unterminated array {name}")


def _doubles(body):
    return [float(tok) for tok in re.findall(r"[-+]?\d+\.\d+e[+-]\d+", body)]


def _ints(body):
    return [int(tok) for tok in re.findall(r"-?\d+", body)]


def load():
    master = _read("boys_coefficients.hpp")
    cheb = _read("boys_coefficients_d_coarse_cheb.hpp")
    rat = _read("boys_coefficients_d_coarse_rat.hpp")

    # kPieces: 66 rows of {a, b, deg, offset}, // Fn
    pieces_text = _array_body(master, "inline constexpr auto kPieces")
    rows = re.findall(
        r"\{\s*([-+0-9.eE]+)\s*,\s*([-+0-9.eE]+)\s*,\s*(\d+)\s*,\s*(\d+)\s*\}", pieces_text
    )
    pieces = [dict(a=float(a), b=float(b), deg=int(d), off=int(o)) for a, b, d, o in rows]

    piece_start = _ints(_array_body(master, "inline constexpr auto kPieceStart"))
    thresholds = _doubles(_array_body(master, "inline constexpr auto kTierThresholds"))

    coeffs = _doubles(_array_body(cheb, "kCoeffs = std::to_array<double>"))

    ext_coeffs = _doubles(_array_body(master, "inline constexpr auto kExtendedBcoeffs"))
    ext_lo = float(re.search(r"inline constexpr double kExtendedBX0 = ([-+0-9.eE]+);", master).group(1))
    x0 = float(re.search(r"inline constexpr double kX0 = ([-+0-9.eE]+);", master).group(1))

    rat_offset = _ints(_array_body(master, "inline constexpr auto kRatAOffset"))
    rat_num_deg = _ints(_array_body(master, "inline constexpr auto kRatANumDeg"))
    rat_den_deg = _ints(_array_body(master, "inline constexpr auto kRatADenDeg"))
    rat_coeffs = _doubles(_array_body(rat, "kRatACoeffs = std::to_array<double>"))

    assert len(pieces) == 66, len(pieces)
    assert len(piece_start) == 34, len(piece_start)
    assert len(coeffs) == 1320, len(coeffs)
    assert len(thresholds) == 33, len(thresholds)
    assert len(ext_coeffs) == 25, len(ext_coeffs)
    assert len(rat_offset) == 66 and len(rat_num_deg) == 66 and len(rat_den_deg) == 66
    assert len(rat_coeffs) == 891, len(rat_coeffs)

    # The offsets must tile the table exactly at deg + 1 each.
    for order in range(33):
        first, last = piece_start[order], piece_start[order + 1]
        expect = first * 0 if False else piece_start[order]
        for i in range(first, last):
            if i + 1 < len(pieces):
                assert pieces[i]["off"] + pieces[i]["deg"] + 1 == pieces[i + 1]["off"], (
                    i,
                    pieces[i],
                    pieces[i + 1],
                )
    assert pieces[-1]["off"] + pieces[-1]["deg"] + 1 == len(coeffs)

    return dict(
        pieces=pieces,
        piece_start=piece_start,
        thresholds=thresholds,
        coeffs=coeffs,
        ext_coeffs=ext_coeffs,
        ext_lo=ext_lo,
        ext_hi=x0,
        rat_offset=rat_offset,
        rat_num_deg=rat_num_deg,
        rat_den_deg=rat_den_deg,
        rat_coeffs=rat_coeffs,
    )


def find_piece(table, order, x):
    first = table["piece_start"][order]
    last = table["piece_start"][order + 1]
    for i in range(first, last - 1):
        if x < table["pieces"][i]["b"]:
            return i
    return last - 1


def coefs_of(table, index):
    piece = table["pieces"][index]
    return table["coeffs"][piece["off"] : piece["off"] + piece["deg"] + 1]


if __name__ == "__main__":
    t = load()
    print("pieces", len(t["pieces"]), "coeffs", len(t["coeffs"]))
    print("bands", t["pieces"][0], t["pieces"][1])
    print("thresholds", t["thresholds"])
    for order in (0, 1, 8, 16, 32):
        i = find_piece(t, order, 0.5)
        print(order, "piece", i, t["pieces"][i], "c0..3", coefs_of(t, i)[:3])
