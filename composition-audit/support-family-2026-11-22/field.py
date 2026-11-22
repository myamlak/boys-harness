"""Exact arithmetic in Q(sqrt 2, sqrt 3) for the mixed-level basis.

Elements are (a,b,c,d) = a + b sqrt2 + c sqrt3 + d sqrt6 with Fraction
coordinates.  The mixed-level orthonormal basis used here lies in this
field.
"""

from __future__ import annotations

from fractions import Fraction

import numpy as np


class F4:
    __slots__ = ("a", "b", "c", "d")

    def __init__(self, a=0, b=0, c=0, d=0):
        self.a = Fraction(a)
        self.b = Fraction(b)
        self.c = Fraction(c)
        self.d = Fraction(d)

    def __add__(self, o):
        if not isinstance(o, F4):
            o = F4(o)
        return F4(self.a + o.a, self.b + o.b, self.c + o.c, self.d + o.d)

    __radd__ = __add__

    def __neg__(self):
        return F4(-self.a, -self.b, -self.c, -self.d)

    def __sub__(self, o):
        if not isinstance(o, F4):
            o = F4(o)
        return F4(self.a - o.a, self.b - o.b, self.c - o.c, self.d - o.d)

    def __rsub__(self, o):
        return F4(o) - self

    def __mul__(self, o):
        if not isinstance(o, F4):
            o = F4(o)
        a, b, c, d = self.a, self.b, self.c, self.d
        e, f, g, h = o.a, o.b, o.c, o.d
        return F4(
            a * e + 2 * b * f + 3 * c * g + 6 * d * h,
            a * f + b * e + 3 * c * h + 3 * d * g,
            a * g + c * e + 2 * b * h + 2 * d * f,
            a * h + d * e + b * g + c * f,
        )

    __rmul__ = __mul__

    def __truediv__(self, o):
        return self * o.inverse()

    def is_zero(self):
        return self.a == 0 and self.b == 0 and self.c == 0 and self.d == 0

    def _mul_matrix(self):
        a, b, c, d = self.a, self.b, self.c, self.d
        # columns: images of 1, sqrt2, sqrt3, sqrt6
        return [
            [a, 2 * b, 3 * c, 6 * d],
            [b, a, 3 * d, 3 * c],
            [c, 2 * d, a, 2 * b],
            [d, c, b, a],
        ]

    def inverse(self):
        if self.is_zero():
            raise ZeroDivisionError
        M = self._mul_matrix()
        n = 4
        # solve M x = e0
        aug = [
            [Fraction(M[i][j]) for j in range(n)] + [Fraction(1 if i == 0 else 0)]
            for i in range(n)
        ]
        for col in range(n):
            piv = next((r for r in range(col, n) if aug[r][col] != 0), None)
            if piv is None:
                raise ZeroDivisionError
            aug[col], aug[piv] = aug[piv], aug[col]
            pv = aug[col][col]
            aug[col] = [v / pv for v in aug[col]]
            for r in range(n):
                if r != col and aug[r][col] != 0:
                    f = aug[r][col]
                    aug[r] = [x - f * y for x, y in zip(aug[r], aug[col])]
        return F4(aug[0][n], aug[1][n], aug[2][n], aug[3][n])

    def __float__(self):
        return (float(self.a) + float(self.b) * 2.0 ** 0.5
                + float(self.c) * 3.0 ** 0.5 + float(self.d) * 6.0 ** 0.5)

    def __repr__(self):
        return (f"F4({self.a},{self.b},{self.c},{self.d})")

    def key(self):
        return (self.a, self.b, self.c, self.d)


ZERO = F4()
ONE = F4(1)


def basis_rows(level_count):
    """Exact rows, same convention as instance.axis_basis."""
    if level_count == 2:
        return [[F4(1), F4(1)], [F4(1), F4(-1)]]
    if level_count == 3:
        # mean (1,1,1); linear sqrt(3/2)(1,0,-1) = (sqrt6/2)(1,0,-1);
        # quadratic (1/sqrt2)(1,-2,1) = (sqrt2/2)(1,-2,1)
        half = Fraction(1, 2)
        return [
            [F4(1), F4(1), F4(1)],
            [F4(0, 0, 0, half), F4(0), F4(0, 0, 0, -half)],
            [F4(0, half), F4(0, -1), F4(0, half)],
        ]
    raise ValueError(level_count)


class Mat:
    """Dense exact matrix over F4."""

    def __init__(self, rows):
        self.rows = [[v if isinstance(v, F4) else F4(v) for v in row] for row in rows]
        self.n = len(self.rows)
        self.m = len(self.rows[0]) if self.rows else 0

    @staticmethod
    def from_array(arr):
        return Mat([[F4(v) for v in row] for row in arr])

    def copy(self):
        return Mat([row[:] for row in self.rows])

    def rank(self):
        A = self.copy()
        r = 0
        n, m = A.n, A.m
        for col in range(m):
            piv = next((i for i in range(r, n) if not A.rows[i][col].is_zero()), None)
            if piv is None:
                continue
            A.rows[r], A.rows[piv] = A.rows[piv], A.rows[r]
            pv = A.rows[r][col]
            A.rows[r] = [v / pv for v in A.rows[r]]
            for i in range(n):
                if i != r and not A.rows[i][col].is_zero():
                    f = A.rows[i][col]
                    A.rows[i] = [x - f * y for x, y in zip(A.rows[i], A.rows[r])]
            r += 1
            if r == n:
                break
        return r

    def inverse(self):
        n = self.n
        if n != self.m:
            raise ValueError("not square")
        aug = [self.rows[i][:] + [ONE if i == j else ZERO for j in range(n)]
               for i in range(n)]
        A = Mat(aug)
        for col in range(n):
            piv = next((i for i in range(col, n) if not A.rows[i][col].is_zero()), None)
            if piv is None:
                raise ZeroDivisionError("singular")
            A.rows[col], A.rows[piv] = A.rows[piv], A.rows[col]
            pv = A.rows[col][col]
            A.rows[col] = [v / pv for v in A.rows[col]]
            for i in range(n):
                if i != col and not A.rows[i][col].is_zero():
                    f = A.rows[i][col]
                    A.rows[i] = [x - f * y for x, y in zip(A.rows[i], A.rows[col])]
        return Mat([row[n:] for row in A.rows])

    def matmul(self, other):
        if self.m != other.n:
            raise ValueError("shape")
        out = []
        for i in range(self.n):
            row = []
            for j in range(other.m):
                s = ZERO
                for k in range(self.m):
                    s = s + self.rows[i][k] * other.rows[k][j]
                row.append(s)
            out.append(row)
        return Mat(out)

    def transpose(self):
        return Mat([[self.rows[i][j] for i in range(self.n)] for j in range(self.m)])

    def submatrix(self, rows, cols):
        return Mat([[self.rows[i][j] for j in cols] for i in rows])

    def scale(self, f):
        return Mat([[v * f for v in row] for row in self.rows])

    def to_float(self):
        return np.array([[float(v) for v in row] for row in self.rows])

    def zero_entries(self):
        return [[v.is_zero() for v in row] for row in self.rows]
