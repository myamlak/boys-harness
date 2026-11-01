"""Decision changes of the frozen-factor surrogates against the true model."""

import sys
from pathlib import Path

from mpmath import mp, mpf

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE / "surrogate"))
import surrogate_control as S  # noqa: E402

mp.dps = 70


def choices(table):
    out = {}
    for key, rows in table.items():
        for row in rows:
            out[(key[0], key[1], row["n"])] = row["state"]
    return out


def main():
    reference, cells, error0, winner, *_ = S.build_measured()
    entries, population, model_from = S.build_cells(reference, cells, error0)
    true_abs = S.predicted(entries, S.A.amplification, "absolute")
    true_rel = S.predicted(entries, S.A.amplification, "relative")
    true1 = S.all_rows(*true_abs, *true_rel, winner, mpf(1))
    true_choice = choices(true1)
    true_set = {k for k, v in true_choice.items() if v != "abstain"}

    def frozen_first_n0(n, x):
        return mpf(1) if n == 0 else mpf(1) / (2 * x)

    def frozen_first(n, x):
        return mpf(1) / (2 * x)

    def unit(n, x):
        return mpf(1)

    for label, fun in (
        ("frozen n0 kept", frozen_first_n0),
        ("frozen n0 also", frozen_first),
        ("unit", unit),
    ):
        f_abs = S.predicted(entries, fun, "absolute")
        f_rel = S.predicted(entries, fun, "relative")
        f1 = S.all_rows(*f_abs, *f_rel, winner, mpf(1))
        fc = choices(f1)
        changed_wrong = []
        changed_abstain = []
        for key in sorted(true_set):
            if fc[key] == "abstain":
                changed_abstain.append(key)
            elif fc[key] != true_choice[key]:
                changed_wrong.append(key)
        changed = changed_wrong + changed_abstain
        print(label, "true commitments", len(true_set),
              "changed", len(changed), "wrong winner", len(changed_wrong),
              "abstentions", len(changed_abstain),
              "all changed relative", all(k[1] == "relative" for k in changed))
        print("  wrong:", changed_wrong)
        print("  abstain:", changed_abstain)


if __name__ == "__main__":
    main()
