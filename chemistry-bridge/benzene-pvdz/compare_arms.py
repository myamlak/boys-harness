"""Compare two bridge arms' dense ERI dumps: max|d|, differing count, energy shift.

Both files are the probe's own text dump: n on the first line, then n^4 lines of
"i j k l value" in the tensor's index order, written at 17 digits (a round-trip).
The comparison streams the value columns in blocks, then, when a matrices file is
given, contracts each tensor with the seed density under the probe's own energy
formula (0.5 sum D D (ij|kl) - 0.25 sum D D (ik|lj)).
"""
import sys
import numpy as np


def stream_values(path, block_bytes=48 * 1024 * 1024):
    with open(path, "rb") as handle:
        first = handle.readline()
        n = int(first)
        carry = b""
        while True:
            chunk = handle.read(block_bytes)
            if not chunk:
                break
            data = carry + chunk
            cut = data.rfind(b"\n")
            if cut < 0:
                carry = data
                continue
            carry = data[cut + 1:]
            block = np.fromstring(data[:cut].decode("ascii"), sep=" ")
            # five numbers per line: i j k l value
            if block.size % 5 != 0:
                raise SystemExit("%s: %d numbers is not a whole number of lines"
                                 % (path, block.size))
            yield block[4::5]
    if carry.strip():
        block = np.fromstring(carry.decode("ascii"), sep=" ")
        yield block[4::5]


def aligned(path_a, path_b, block=1 << 22):
    """Yields equal-length value blocks from both files.

    The files' lines differ in width because their values do, so a fixed byte
    block splits into a different number of values on each side; the buffers
    below re-align the two streams before anything is compared.
    """
    it_a, it_b = iter(stream_values(path_a)), iter(stream_values(path_b))
    buf_a = np.empty(0)
    buf_b = np.empty(0)
    done_a = done_b = False
    while True:
        while buf_a.size < block and not done_a:
            try:
                buf_a = np.concatenate([buf_a, next(it_a)])
            except StopIteration:
                done_a = True
        while buf_b.size < block and not done_b:
            try:
                buf_b = np.concatenate([buf_b, next(it_b)])
            except StopIteration:
                done_b = True
        if buf_a.size == 0 and buf_b.size == 0:
            return
        take = min(buf_a.size, buf_b.size)
        if take == 0:
            raise SystemExit("the two files hold a different number of values")
        yield buf_a[:take], buf_b[:take]
        buf_a = buf_a[take:]
        buf_b = buf_b[take:]
        if buf_a.size == 0 and buf_b.size == 0 and done_a and done_b:
            return


def compare(path_a, path_b):
    max_abs = 0.0
    differing = 0
    total = 0
    for a, b in aligned(path_a, path_b):
        d = np.abs(a - b)
        max_abs = max(max_abs, float(d.max()))
        differing += int(np.count_nonzero(d))
        total += a.size
    return max_abs, differing, total


def read_tensor(path):
    with open(path, "rb") as handle:
        n = int(handle.readline())
    values = np.concatenate(list(stream_values(path)))
    return n, values.reshape((n, n, n, n))


def read_density(matrices_path, name):
    with open(matrices_path) as handle:
        lines = handle.read().split("\n")
    start = None
    for index, line in enumerate(lines):
        if line.startswith("matrix " + name + " "):
            start = index
            break
    if start is None:
        raise SystemExit("matrix %s not found in %s" % (name, matrices_path))
    _, _, rows, cols = lines[start].split()
    rows, cols = int(rows), int(cols)
    block = np.array([float(token) for token in
                      " ".join(lines[start + 1:start + 1 + rows]).split()])
    return block.reshape((rows, cols))


def two_electron(eri, density):
    n = eri.shape[0]
    dvec = density.reshape(-1)
    coulomb = 0.5 * float(dvec @ (eri.reshape((n * n, n * n)) @ dvec))
    exchanged = np.ascontiguousarray(eri.transpose(0, 2, 3, 1))
    exchange = -0.25 * float(dvec @ (exchanged.reshape((n * n, n * n)) @ dvec))
    return coulomb + exchange


def main():
    if len(sys.argv) < 3:
        raise SystemExit("usage: compare_arms.py <arm-a.eri.txt> <arm-b.eri.txt> [matrices.txt <density>]")
    max_abs, differing, total = compare(sys.argv[1], sys.argv[2])
    print("values %d | max|d| %.17e | differing %d" % (total, max_abs, differing))
    if len(sys.argv) >= 5:
        n_a, eri_a = read_tensor(sys.argv[1])
        n_b, eri_b = read_tensor(sys.argv[2])
        if n_a != n_b:
            raise SystemExit("tensor sizes differ")
        density = read_density(sys.argv[3], sys.argv[4])
        e_a = two_electron(eri_a, density)
        e_b = two_electron(eri_b, density)
        print("two-electron energy  a %.17e | b %.17e | shift %.3e"
              % (e_a, e_b, e_b - e_a))


if __name__ == "__main__":
    main()
