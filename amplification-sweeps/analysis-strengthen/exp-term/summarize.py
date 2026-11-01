#!/usr/bin/env python3
"""Every count the model's section states, read from the regenerated records.

Inputs (read only):
    analysis-strengthen/out/{second,first}/prediction.csv   the new model
    analysis-strengthen/out/{second,first}/scan_all.csv     the ladder
    analysis-strengthen/exp-term/out/v5.2.0-exp0/...        the previous model

Prints the per-band and aggregate commitment counts, the confusion counts,
the coverage over the 100 comparisons that are not same-fit abstentions, the
always-direct score over the whole set and over the model's committed set,
the ladder on both sweeps, and the row-by-row changes from the previous
model.
"""

import csv
from pathlib import Path

HERE = Path(__file__).resolve().parent
STRENGTHEN = HERE.parent
OUT = STRENGTHEN / "out"
OLD = HERE / "out" / "v5.2.0-exp0"

BANDS = ("A1", "A2")
METRICS = ("absolute", "relative")
NMAX = 32
SAME_FIT = {("A1", "absolute"), ("A1", "relative")}  # orders 17..32


def read_prediction(path):
    rows = {}
    with open(path, newline="") as handle:
        for row in csv.DictReader(handle):
            rows[(row["band"], row["metric"], int(row["order"]))] = row
    return rows


def committed(row):
    return row["state"] != "abstain"


def main():
    for sweep in ("second", "first"):
        new = read_prediction(OUT / sweep / "prediction.csv")
        old = read_prediction(OLD / sweep / "prediction.csv")
        print("=== %s sweep" % sweep)
        total = [0, 0, 0]
        for band in BANDS:
            for metric in METRICS:
                keys = [(band, metric, n) for n in range(NMAX + 1)]
                commits = sum(1 for k in keys if committed(new[k]))
                right = sum(1 for k in keys if committed(new[k]) and new[k]["right_where_committed"] == "1")
                abstains = len(keys) - commits
                rec = sum(1 for k in keys if committed(new[k]) and new[k]["choice"] == "ship")
                dir_ = sum(1 for k in keys if committed(new[k]) and new[k]["choice"] == "product")
                won_rec = sum(1 for k in keys if committed(new[k]) and new[k]["choice"] == "ship" and new[k]["right_where_committed"] == "1")
                won_dir = sum(1 for k in keys if committed(new[k]) and new[k]["choice"] == "product" and new[k]["right_where_committed"] == "1")
                lost = commits - right
                direct_wins = sum(1 for k in keys if new[k]["measured_winner"] == "product")
                ties = sum(1 for k in keys if new[k]["measured_winner"] == "tie")
                print("  %s %-8s commits %2d right %2d abstains %2d | rec %2d (won %2d) dir %2d (won %2d) lost %d | direct wins %2d ties %d"
                      % (band, metric, commits, right, abstains, rec, won_rec, dir_, won_dir, lost, direct_wins, ties))
                total[0] += commits
                total[1] += right
                total[2] += abstains
        print("  all: commits %d right %d abstains %d" % tuple(total))

        # coverage over the 100 comparisons outside the same-fit block
        rest_keys = []
        same_keys = []
        for band in BANDS:
            for metric in METRICS:
                for n in range(NMAX + 1):
                    if (band, metric) in SAME_FIT and n >= 17:
                        same_keys.append((band, metric, n))
                    else:
                        rest_keys.append((band, metric, n))
        rest_commits = sum(1 for k in rest_keys if committed(new[k]))
        rest_right = sum(1 for k in rest_keys if committed(new[k]) and new[k]["right_where_committed"] == "1")
        rest_dir_right = sum(1 for k in rest_keys if new[k]["measured_winner"] == "product")
        rest_dir_wrong = sum(1 for k in rest_keys if new[k]["measured_winner"] == "ship")
        rest_ties = sum(1 for k in rest_keys if new[k]["measured_winner"] == "tie")
        print("  same-fit block: %d keys, commits %d"
              % (len(same_keys), sum(1 for k in same_keys if committed(new[k]))))
        print("  remaining 100: commits %d right %d; always-direct right %d wrong %d tied %d"
              % (rest_commits, rest_right, rest_dir_right, rest_dir_wrong, rest_ties))
        all_dir_right = sum(1 for k in new if new[k]["measured_winner"] == "product")
        all_ties = sum(1 for k in new if new[k]["measured_winner"] == "tie")
        print("  all 132: always-direct right %d, tied %d" % (all_dir_right, all_ties))

        # always-direct on the model's committed set
        dir_on_commits = sum(1 for k in new if committed(new[k]) and new[k]["measured_winner"] == "product")
        print("  always-direct right on the model's committed set: %d" % dir_on_commits)

        changed = []
        for key in sorted(new):
            a = (old[key]["state"], old[key]["choice"])
            b = (new[key]["state"], new[key]["choice"])
            if a != b:
                changed.append((key, a, b))
        print("  changed rows from the previous model:")
        for key, a, b in changed:
            print("    %s %-8s n=%2d %s/%s -> %s/%s"
                  % (key[0], key[1], key[2], a[0], a[1] or "-", b[0], b[1] or "-"))
        print()

    print("=== ladder (aggregate commits/right/wrong/abstains)")
    for sweep in ("second", "first"):
        print("  %s:" % sweep)
        with open(OUT / sweep / "scan_all.csv", newline="") as handle:
            for row in csv.DictReader(handle):
                print("    margin %-4s %3s commits %3s right %3s wrong %3s abstains"
                      % (row["margin"], row["commits"], row["right"], row["wrong"], row["abstains"]))


if __name__ == "__main__":
    main()
