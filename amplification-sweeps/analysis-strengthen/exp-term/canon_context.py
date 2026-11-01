"""Print canonical document lines around a search term, for exact ledger quotes."""

import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
PAPER = HERE.parents[1]
sys.path.insert(0, str(PAPER / "tools"))
import check_paper_numbers as C  # noqa: E402

DOCS = {
    "main.tex": C.Document(PAPER / "main.tex", True),
    "supplement/main.tex": C.Document(PAPER / "supplement" / "main.tex", True),
    "cover-letter.md": C.Document(PAPER / "cover-letter.md", False),
}


def main():
    doc_name, term = sys.argv[1], sys.argv[2]
    before = int(sys.argv[3]) if len(sys.argv) > 3 else 1
    after = int(sys.argv[4]) if len(sys.argv) > 4 else 1
    canon = DOCS[doc_name].text.split("\n")
    hits = [i for i, line in enumerate(canon) if re.search(term, line, re.IGNORECASE)]
    shown = set()
    for i in hits:
        for j in range(max(0, i - before), min(len(canon), i + after + 1)):
            if j in shown:
                continue
            shown.add(j)
            print("%4d| %s" % (j, canon[j]))
        print("    " + "-" * 60)
    if not hits:
        print("NO MATCH for", term)


if __name__ == "__main__":
    main()
