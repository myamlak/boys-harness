#!/usr/bin/env python3
"""Synthetic self-test of pair_analyze.py (scratch copy only)."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import pair_analyze as pa

OUT = Path(r"C:/Users/myaml/dev/bladawiec-pair/build-pair/selftest")

states = pa.model_states()
arguments = [(2, 1, 0.5), (2, 1, 6.0), (1, 1, 12.5), (0, 1, 0.0)]

OUT.parent.mkdir(parents=True, exist_ok=True)

for arm in pa.ARMS:
    with open(OUT.with_suffix(".%s.args.txt" % arm), "w") as handle:
        handle.write("# arm=%s nmax count x paper-a direct engine\n" % arm)
        for nmax, count, x in arguments:
            paper_a = [1.0 / (2 * n + 1) * (1.0 - 0.1 * x) for n in range(nmax + 1)]
            direct = [1.0 / (2 * n + 1) * (1.0 - 0.05 * x) for n in range(nmax + 1)]
            engine = []
            for order in range(nmax + 1):
                column = pa.selected_column(arm, x, order, states)
                engine.append(paper_a[order] if column == "paper-a" else direct[order])
            tokens = [str(nmax), str(count), repr(x)]
            tokens += [repr(v) for v in paper_a]
            tokens += [repr(v) for v in direct]
            tokens += [repr(v) for v in engine]
            handle.write(" ".join(tokens) + "\n")
    with open(OUT.with_suffix(".%s.summary.txt" % arm), "w") as handle:
        handle.write("seed_e2e=1.0\ndefault_e2e=1.0\ntight_e2e=1.0\n")
        handle.write("scf_default_total_energy=-74.0\nscf_tight_total_energy=-74.0\n")
    with open(OUT.with_suffix(".%s.eri.txt" % arm), "w") as handle:
        handle.write("1\n0 0 0 0 %r\n" % (1.0 + pa.ARMS.index(arm) * 1e-16))

result = OUT.with_suffix(".analysis.txt")
sys.argv = ["pair_analyze.py", str(OUT), str(result)]
pa.main()

text = result.read_text()
for needle in ("referee's bar", "leave-one-out", "tensors"):
    if needle not in text:
        raise SystemExit("missing section: " + needle)
print("SELFTEST OK")
