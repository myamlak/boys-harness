"""Compare the platform exponential with a 70-digit reference.

Input: exp/exp_bits.txt (xbits expbits halfbits).  Prints the worst relative
error in units of 2^-53, the count above 0.5 and 1, and whether Python's
math.exp agrees bit for bit with the platform probe.
"""
import math
import struct
from mpmath import mp, mpf

mp.dps = 70
U = mpf(2) ** -53
HERE = r"C:/Users/myaml/dev/_a-cert-run/exp"

worst = (mpf(0), None)
rel_worst = (mpf(0), None)
over_half = over_one = 0
py_mismatch = 0
half_mismatch = 0
count = 0
for line in open(HERE + "/exp_bits.txt"):
    xb_s, eb_s, hb_s = line.split()
    xb, eb, hb = int(xb_s, 16), int(eb_s, 16), int(hb_s, 16)
    x = struct.unpack(">d", struct.pack(">Q", xb))[0]
    e = struct.unpack(">d", struct.pack(">Q", eb))[0]
    h = struct.unpack(">d", struct.pack(">Q", hb))[0]
    ref = mp.exp(-mpf(x))
    ulps = abs(mpf(e) - ref) / ref / U
    rel = abs(mpf(e) - ref) / ref
    count += 1
    if ulps > worst[0]:
        worst = (ulps, (x, e))
    if rel > rel_worst[0]:
        rel_worst = (rel, (x, e))
    if ulps > mpf("0.5"):
        over_half += 1
    if ulps > 1:
        over_one += 1
    if math.exp(-x) != e:
        py_mismatch += 1
    if h != e * 0.5:
        half_mismatch += 1

print("arguments", count)
print("worst |delta|/ref in 2^-53 units:", mp.nstr(worst[0], 12), "at x =", worst[1][0])
print("worst relative:", mp.nstr(rel_worst[0], 12))
print("arguments above 0.5 unit:", over_half, "above 1 unit:", over_one)
print("python math.exp mismatches:", py_mismatch)
print("half-scaling mismatches:", half_mismatch)
print("conditional eta_E <= 2^-52 supported:", worst[0] < 2)
