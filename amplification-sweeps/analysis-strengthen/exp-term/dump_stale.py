import json
from pathlib import Path

root = Path(__file__).resolve().parents[2]
led = json.loads((root / "tools" / "ledger.json").read_text(encoding="utf-8"))
claims = {c["id"]: c for c in led["claims"]}
ids = """abstract-band-ratios abstract-second-size abstract-commits abstract-abstains
abstract-first-commits abstract-first-right abstract-oracle intro-band-ratios
intro-commits intro-abstains intro-first-commits intro-first-right intro-oracle
intro-bases intro-cell-win intro-sto-margin prediction-b1-rel prediction-b2-rel
prediction-total replication-commits replication-pop-right replication-whole-right
replication-direct replication-ladder samefit-coverage group-wins group-commits
below-theta0 surrogate-replay surrogate-frozen-a1 surrogate-an-one surrogate-rise
b1-abstentions frontier-half frontier-term frontier-tail interval-cp interval-lost
interval-shape interval-first e2e-fallback conclusion-e2e bound-full
s8-commits-right s8-lost s8-retuned-margin s8-ceiling-first-wrong s8-true-range
s8-six-decisions s8-a1 s8-an s8-rise-supplies s8-a1-breakdown s8-an-breakdown
s12-first-commits s12-first-population-right s12-first-whole-right s12-second-scores
s12-direct-right cover-commits cover-abstains cover-lost cover-apriori-certifies
abstract-reference-precision abstract-abstains-fallback intro-abstains-fallback
intro-report-arguments results-truncation-16 results-remaining-51 results-charge-92
results-67-of-67 results-caption-totals results-model-score limitations-certifies
cover-reference-precision main-first-right-68 s12-67-of-67""".split()
out = []
for i in ids:
    c = claims[i]
    out.append("=" * 8 + " " + i + " | " + str(c.get("file")) + " | " + str(c.get("location")))
    out.append("QUOTE: " + repr(c.get("quote")))
    out.append("TEMPLATE: " + repr(c.get("template")))
    out.append("VALUES: " + json.dumps(c.get("values")))
print("\n".join(out))
