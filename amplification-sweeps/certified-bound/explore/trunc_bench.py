"""Timing and width check of the certified truncation interval."""
import argparse
import sys
import time
from pathlib import Path

from mpmath import mp, mpf

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "code"))
import record_paths

_parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
record_paths.add_inputs(_parser, deposit_arg=True)
_args = _parser.parse_args()
HERE = record_paths.deposit(_args.deposit)
with record_paths.deposit_on_path(HERE):
    import bounds
    import exact

mp.dps = 60

x = 6.3
for n, bi in ((5, 1), (20, 0)):
    e = exact.fit_exact_error(n, bi, x)
    print("n", n, "exact", mp.nstr(e, 10))
    for samples in (16, 32, 64, 128, 400, 1600):
        t0 = time.time()
        lo, hi = bounds.fit_error(n, bi, x, "sampled", samples)
        dt = time.time() - t0
        print("  samples %4d  lo %s  hi %s  width %s  /exact %s  %.3fs"
              % (samples, mp.nstr(lo, 8), mp.nstr(hi, 8),
                 mp.nstr(hi - lo, 8),
                 mp.nstr((hi - lo) / abs(e), 6) if e else "-", dt))
t0 = time.time()
for _ in range(10):
    bounds.seed_error(mpf(6.3), "sampled", 32)
print("seed_error x10 samples32 %.3fs" % (time.time() - t0))
