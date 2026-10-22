"""Send query batches to the crlibm driver running under WSL.

The driver (work/driver.c, built against crlibm 1.0beta5 as configured for
x86-64 Linux) reads "func mode hexin" lines on stdin and writes
"func mode hexin hexout" lines on stdout.
"""
import subprocess
import struct

WSL = ["wsl.exe", "-e", "bash", "-lc", "~/crdrv"]


def double_to_hex(d):
    return struct.pack(">d", d).hex()


def hex_to_double(h):
    return struct.unpack(">d", bytes.fromhex(h))[0]


def evaluate(queries, batch=200000):
    """queries: iterable of (func, mode, inhex) or (func, 'rn', hex, hex2).

    Returns a list of output hex strings (or 'ERR'/'nan' markers) in order.
    """
    queries = list(queries)
    out = []
    for i in range(0, len(queries), batch):
        chunk = queries[i:i + batch]
        lines = []
        for q in chunk:
            if len(q) == 4:
                lines.append("%s %s %s %s" % (q[0], q[1], q[2], q[3]))
            else:
                lines.append("%s %s %s" % (q[0], q[1], q[2]))
        p = subprocess.run(WSL, input="\n".join(lines) + "\n",
                           capture_output=True, text=True)
        res = p.stdout.split("\n")
        seen = 0
        for line in res:
            f = line.split()
            if not f:
                continue
            if len(f) >= 4:
                out.append(f[-1])
                seen += 1
        if seen != len(chunk):
            raise RuntimeError("driver returned %d lines for %d queries"
                               % (seen, len(chunk)))
    return out
