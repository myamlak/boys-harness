"""Timing and width check of the certified truncation interval."""
import sys
import time

sys.path.insert(0, r"C:/Users/myaml/dev/_harness-v1023/amplification-sweeps/apriori-bound")
from mpmath import mp, mpf

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
