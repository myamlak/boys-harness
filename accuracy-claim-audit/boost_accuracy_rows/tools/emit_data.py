"""Write the published test points of a data group out as the gate's sampling data.

Each Boost.Math data file holds one array per named group.  The layout differs by family:
erf and gamma files carry (x, f(x)), while the Bessel files carry (order, x, value).  The gate
reads x and, where a row's order is taken from the data, the order; so this normalises both
layouts to "x[, order]".
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from survey import arrays_in  # noqa: E402

# file, array, where the x column is, where the order column is (None if the group fixes it)
PLAN = [
    ("erf_small_data.ipp", "erf_small_data", 0, None),
    ("erf_data.ipp", "erf_data", 0, None),
    ("erf_large_data.ipp", "erf_large_data", 0, None),
    ("expinti_data.ipp", "expinti_data", 0, None),
    ("digamma_neg_data.ipp", "digamma_neg_data", 0, None),
    ("digamma_data.ipp", "digamma_data", 0, None),
    ("test_gamma_data.ipp", "near_1", 0, None),
    ("test_gamma_data.ipp", "near_m10", 0, None),
    ("bessel_j_data.ipp", "bessel_j_data", 1, 0),
    ("bessel_i_data.ipp", "bessel_i_data", 1, 0),
]


def main() -> int:
    test_dir = Path(sys.argv[1])
    out_dir = Path(sys.argv[2])
    out_dir.mkdir(parents=True, exist_ok=True)
    for file_name, array, x_col, nu_col in PLAN:
        found = False
        for name, _, grid in arrays_in(test_dir / file_name):
            if name != array:
                continue
            found = True
            dest = out_dir / f"{array}.csv"
            with dest.open("w", encoding="utf-8") as fh:
                for row in grid:
                    x = row[x_col]
                    if nu_col is None:
                        fh.write(f"{x!r}\n")
                    else:
                        fh.write(f"{x!r},{row[nu_col]!r}\n")
            print(f"{dest.name:22s} {len(grid):5d} points from {file_name}")
        if not found:
            print(f"  MISSING: {array} not found in {file_name}")
            return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
