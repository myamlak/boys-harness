#!/usr/bin/env python3
"""Exact algebraic check of the fused-fp32 component-120 share at n=2, x=4.

The orthonormal contrast entries lie in Q(sqrt2, sqrt3).  Every recorded
F0-scaled ratio is a binary rational.  Each coefficient, the non-mean total and
the share are therefore exact elements of Q(sqrt2, sqrt3, sqrt6), which this
script carries as four Fractions.  No floating point is used in the decision.

Checks whether 10 c_120^2 - D0 is exactly zero, i.e. whether the share is
exactly one tenth.
"""
from __future__ import annotations

import importlib.util
import sys
from fractions import Fraction
from pathlib import Path

import numpy as np

SCRIPT = Path(__file__).resolve().parent / "fullcube_decomposition.py"
FUSED = (Path(__file__).resolve().parents[1]
         / "records" / "cells" / "cells-msvc-fused.csv.xz")
SEPARATE = (Path(__file__).resolve().parents[1]
            / "records" / "cells" / "cells-msvc-separate.csv.xz")


def load_module(path: Path):
    spec = importlib.util.spec_from_file_location("fcd", str(path))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class F:
    """a + b*sqrt2 + c*sqrt3 + d*sqrt6 with rational components."""

    __slots__ = ("a", "b", "c", "d")

    def __init__(self, a=0, b=0, c=0, d=0):
        self.a = Fraction(a)
        self.b = Fraction(b)
        self.c = Fraction(c)
        self.d = Fraction(d)

    def __add__(self, other):
        other = other if isinstance(other, F) else F(other)
        return F(self.a + other.a, self.b + other.b, self.c + other.c, self.d + other.d)

    __radd__ = __add__

    def __neg__(self):
        return F(-self.a, -self.b, -self.c, -self.d)

    def __sub__(self, other):
        return self + (-other if isinstance(other, F) else F(-other))

    def __mul__(self, other):
        other = other if isinstance(other, F) else F(other)
        a, b, c, d = self.a, self.b, self.c, self.d
        e, f, g, h = other.a, other.b, other.c, other.d
        return F(
            a * e + 2 * b * f + 3 * c * g + 6 * d * h,
            a * f + b * e + 3 * c * h + 3 * d * g,
            a * g + c * e + 2 * b * h + 2 * d * f,
            a * h + d * e + b * g + c * f,
        )

    __rmul__ = __mul__

    def is_zero(self):
        return self.a == 0 and self.b == 0 and self.c == 0 and self.d == 0

    def __repr__(self):
        return f"F({self.a},{self.b},{self.c},{self.d})"


def phi_matrix(fcd):
    axes = []
    for count in fcd.LEVEL_COUNTS:
        if count == 2:
            axes.append([[F(1), F(1)], [F(1), F(-1)]])
        else:
            s6_half = F(0, 0, 0, Fraction(1, 2))      # sqrt(6)/2 = sqrt(3/2)
            s2_half = F(0, Fraction(1, 2), 0, 0)      # sqrt(2)/2 = 1/sqrt(2)
            s2_one = F(0, 1, 0, 0)                    # sqrt(2) = 2/sqrt(2)
            axes.append(
                [
                    [F(1), F(1), F(1)],
                    [s6_half, F(0), -s6_half],
                    [s2_half, -s2_one, s2_half],
                ]
            )
    grid = axes[0]
    for matrix in axes[1:]:
        grid = [[x * y for x in row_x for y in row_y] for row_x in grid for row_y in matrix]
    return grid  # [component][cell]


def exact_coefficients(row, phi):
    """c_k = (1/144) sum_i r_i phi_k,i, exact."""
    r = [F(Fraction(float(v))) for v in row]
    out = []
    for component in range(len(phi)):
        total = F(0)
        row_phi = phi[component]
        for cell in range(len(r)):
            total = total + r[cell] * row_phi[cell]
        out.append(F(total.a / 144, total.b / 144, total.c / 144, total.d / 144))
    return out


def main() -> int:
    fcd = load_module(SCRIPT)
    phi = phi_matrix(fcd)
    print(f"exact Phi matrix: {len(phi)} x {len(phi[0])}")
    cases = [
        ("fused", FUSED, "fp32", 2, 4.0, 120),
        ("separate", SEPARATE, "fp32", 0, 0.25, 120),
    ]
    for route, path, lane, n, x, component in cases:
        run = fcd.load_route(str(path), route)
        lane_index = fcd.LANES.index(lane)
        point = None
        for p in range(run["R"].shape[1]):
            loc = fcd.point_location(run, p)
            if loc["n"] == n and loc["x"] == x:
                point = p
                break
        row = run["R"][lane_index, point]
        c = exact_coefficients(row, phi)
        d0 = F(0)
        for k in range(1, len(c)):
            d0 = d0 + c[k] * c[k]
        ten_c2_minus_d0 = (c[component] * c[component]) * 10 - d0
        print()
        print(f"{route} {lane} n={n} x={x} component {component}")
        print(f"  unique recorded ratios: {np.unique(row).size}")
        print(f"  c_{component} exact field: {c[component]}")
        print(f"  D0 exact field: {d0}")
        print(f"  10*c^2 - D0 exact field: {ten_c2_minus_d0}")
        print(f"  share exactly 1/10: {ten_c2_minus_d0.is_zero()}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
