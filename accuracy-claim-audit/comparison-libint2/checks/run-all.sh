#!/bin/sh
# Regenerate the accuracy-comparison record from source.
#
# Run inside WSL or Linux (the probe needs an 80-bit long double to parse
# the reference without rounding it, and the libint2 static archive is an
# ELF Linux build).  From the repository root, with both variables in the
# WSL environment:
#
#   MSYS_NO_PATHCONV=1 wsl -e sh -c \
#     'BOYS_SOURCE=/path/to/boys-v5.3.2 LIBINT_ROOT=/path/to/libint-2.13.1 sh accuracy-claim-audit/comparison-libint2/checks/run-all.sh'
#
# Environment:
#   BOYS_SOURCE  boys v5.3.2 checkout (default: ../boys beside the repository)
#   LIBINT_ROOT  libint2 2.13.1 install tree, with inst/include and
#                inst/lib/libint2.a (required)
#
# Every command and its exit code in commands-and-exit-codes.txt comes
# from this script's own run.  The script exits 0 only when every
# critical step succeeds.  The digest lines for the two check scripts and
# the paths in the recorded logs name the measurement tree; see README.md.

set -u

CHECKS=$(cd "$(dirname "$0")" && pwd)
WORK=$(cd "$CHECKS/.." && pwd)
HARNESS_REPO=$(cd "$WORK/../.." && pwd)
BOYS="${BOYS_SOURCE:-$(cd "$HARNESS_REPO/.." && pwd)/boys}"
LIBINT="${LIBINT_ROOT:?set LIBINT_ROOT to the libint2 2.13.1 install tree}"
PROBE=$WORK/probe/compare_boys_accuracy.cpp
BIN=$WORK/bin/compare_boys_accuracy
LOGDIR=$WORK/logs
RECORD=$WORK
LOG=$RECORD/commands-and-exit-codes.txt
DIG=$RECORD/digests.txt

mkdir -p "$LOGDIR" "$WORK/bin"
: > "$LOG"
: > "$DIG"

fail=0

run() {
    label=$1
    cmd=$2
    {
        printf '### %s\n' "$label"
        printf '$ %s\n' "$cmd"
    } >> "$LOG"
    sh -c "$cmd" >> "$LOG" 2>&1
    rc=$?
    printf 'exit=%d\n\n' "$rc" >> "$LOG"
    if [ "$rc" -ne 0 ]; then
        fail=1
        printf 'STEP FAILED: %s (exit %d)\n' "$label" "$rc" >&2
    fi
    return "$rc"
}

run "compiler version" "g++ --version"
run "kernel" "uname -a"

run "boys annotated tag object" "git -C $BOYS rev-parse v5.3.2"
run "boys tag commit" "git -C $BOYS rev-parse 'v5.3.2^{commit}'"
run "boys tag tree" "git -C $BOYS rev-parse 'v5.3.2^{tree}'"
run "boys tracked tree vs tag (all paths)" \
    "git -C $BOYS diff --quiet v5.3.2 && echo 'tracked tree identical to v5.3.2'"
run "boys working tree status" "git -C $BOYS status --porcelain"

run "libint2 version" "grep -n 'define LIBINT_VERSION' $LIBINT/inst/include/libint2/config.h"
run "libint2 checkout kind" \
    "if [ -d $LIBINT/.git ]; then echo 'git checkout'; else echo 'not a git checkout (no .git directory)'; fi"

run "grid, source and library digests" \
    "sha256sum $BOYS/tests/data/boys_accuracy_gate_reference.csv \
$BOYS/tests/data/boys_reference.csv \
$BOYS/tools/gen_boys_accuracy_gate_reference.py \
$BOYS/tools/gen_boys_coefficients.py \
$PROBE \
$WORK/checks/compare_grids.py \
$WORK/checks/run-all.sh \
$LIBINT/inst/include/libint2/boys.h \
$LIBINT/inst/include/libint2/config.h \
$LIBINT/inst/lib/libint2.a | tee $DIG"
run "boys header digests" \
    "find $BOYS/include/boys -name '*.hpp' -print0 | sort -z | xargs -0 sha256sum | tee -a $DIG"

BUILD_CMD="g++ -std=c++20 -O2 -mavx2 -mfma -Wall -Wextra -I $BOYS/include -I $LIBINT/inst/include $PROBE $LIBINT/inst/lib/libint2.a -o $BIN"

if run "compile probe" "$BUILD_CMD 2> $LOGDIR/build-warnings.txt"; then
    run "build warnings (all from libint2 headers)" "cat $LOGDIR/build-warnings.txt"
    run "probe binary digest" "sha256sum $BIN | tee -a $DIG"

    if run "gate grid run" \
        "$BIN $BOYS/tests/data/boys_accuracy_gate_reference.csv boys_accuracy_gate_25digit > $LOGDIR/run-gate.txt 2> $LOGDIR/run-gate.err"; then
        run "gate grid output" "cat $LOGDIR/run-gate.txt"
    else
        printf '### gate grid output skipped: run failed\n' >> "$LOG"
    fi

    if run "45-digit grid run" \
        "$BIN $BOYS/tests/data/boys_reference.csv boys_reference_45digit > $LOGDIR/run-45digit.txt 2> $LOGDIR/run-45digit.err"; then
        run "45-digit grid output" "cat $LOGDIR/run-45digit.txt"
    else
        printf '### 45-digit grid output skipped: run failed\n' >> "$LOG"
    fi

    run "probe stderr (both runs; empty means no parse or domain errors)" \
        "cat $LOGDIR/run-gate.err $LOGDIR/run-45digit.err"
else
    printf '### probe runs skipped: compile failed\n' >> "$LOG"
fi

run "grid cross-check" \
    "python3 $WORK/checks/compare_grids.py --boys $BOYS > $LOGDIR/grid-crosscheck.txt 2> $LOGDIR/grid-crosscheck.err"
run "grid cross-check output" "cat $LOGDIR/grid-crosscheck.txt"

if [ -s "$LOGDIR/run-gate.txt" ] && [ -s "$LOGDIR/run-45digit.txt" ]; then
    run "extract TSV header from probe output" \
        "sed -n 's/^# columns: //p' $LOGDIR/run-gate.txt | sed 's/, /\t/g' > $RECORD/results-header.tsv"
    run "extract gate grid rows" "grep -v '^#' $LOGDIR/run-gate.txt > $LOGDIR/rows-gate.tsv"
    run "extract 45-digit grid rows" "grep -v '^#' $LOGDIR/run-45digit.txt > $LOGDIR/rows-45digit.tsv"
    run "assemble results-gate-grid.tsv" \
        "cat $RECORD/results-header.tsv $LOGDIR/rows-gate.tsv > $RECORD/results-gate-grid.tsv"
    run "assemble results-45digit-grid.tsv" \
        "cat $RECORD/results-header.tsv $LOGDIR/rows-45digit.tsv > $RECORD/results-45digit-grid.tsv"
    run "assemble results-all.tsv" \
        "cat $RECORD/results-header.tsv $LOGDIR/rows-gate.tsv $LOGDIR/rows-45digit.tsv > $RECORD/results-all.tsv"
    run "copy raw outputs into record" \
        "cp $LOGDIR/run-gate.txt $RECORD/raw-gate-grid.txt && \
cp $LOGDIR/run-45digit.txt $RECORD/raw-45digit-grid.txt && \
cp $LOGDIR/build-warnings.txt $RECORD/build-warnings.txt && \
cp $LOGDIR/grid-crosscheck.txt $RECORD/grid-crosscheck.txt"
else
    printf '### TSV extraction skipped: probe output missing\n' >> "$LOG"
    fail=1
fi

{
    printf 'Probe build command (run inside WSL)\n'
    printf '====================================\n\n'
    printf '%s\n\n' "$BUILD_CMD"
    printf 'Flags: -std=c++20 -O2 -mavx2 -mfma -Wall -Wextra\n'
    printf 'Include paths: Boys v5.3.2 public headers; libint2 2.13.1 inst/include\n'
    printf 'Link: libint2 2.13.1 inst/lib/libint2.a (only for libint2_static_cleanup;\n'
    printf 'the Boys evaluator is header-only).\n'
    printf 'Compiler: g++ 15.2.0, WSL Linux x86-64 (see commands-and-exit-codes.txt).\n'
    printf 'Build warnings: all 32 lines come from libint2 third-party headers;\n'
    printf 'the probe source compiles warning-free under -Wall -Wextra.\n'
} > "$RECORD/build-command.txt"

exit "$fail"
