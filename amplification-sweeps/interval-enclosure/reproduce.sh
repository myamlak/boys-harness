#!/bin/sh
# Reproduce the Q-417 interval/enclosure record from the deposit.
#
# Run from anywhere; the script changes to its own directory.  The input
# records are read from the repository root two levels above that directory,
# or from BOYS_HARNESS when set:
#
#   sh reproduce.sh
#   BOYS_HARNESS=/path/to/boys-harness sh reproduce.sh
#
# Requirements: Python 3.10 or newer with mpmath.  No checkout of the boys
# library is read; the deposited input records carry their own provenance.
# The eight commands and their measured exit codes are in logs/commands.log.

set -eu
cd "$(dirname "$0")"
export PYTHONDONTWRITEBYTECODE=1

python -B work/q417_interval.py --case first9353
python -B work/q417_interval.py --case first9352
python -B work/q417_interval.py --case second
python -B work/q417_interval.py --case heldout
python -B work/q417_containment.py --case second
python -B work/q417_containment.py --case heldout
python -B work/q417_containment.py --case first9353
python -B work/q417_checks.py
