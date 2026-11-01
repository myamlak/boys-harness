import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
PAPER = HERE.parents[1]
sys.path.insert(0, str(PAPER / "tools"))
import check_paper_numbers as C  # noqa: E402

led = json.loads((PAPER / "tools" / "ledger.json").read_text(encoding="utf-8"))
docs = {
    "main.tex": C.Document(PAPER / "main.tex", True),
    "supplement/main.tex": C.Document(PAPER / "supplement" / "main.tex", True),
    "cover-letter.md": C.Document(PAPER / "cover-letter.md", False),
}
for c in led["claims"]:
    rel = c.get("file", "main.tex")
    doc = docs[rel]
    matches = list(C.quote_regex(c.get("quote", "")).finditer(doc.text))
    occ = int(c.get("occurrence", 0))
    status = "OK" if len(matches) > occ else ("SHIFT" if matches else "GONE")
    print("%-28s %-20s occ=%d matches=%d %s | %s" % (
        c["id"], rel, occ, len(matches), status, c.get("quote", "")[:70]))
