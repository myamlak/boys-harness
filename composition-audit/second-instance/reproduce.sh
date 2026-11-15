#!/usr/bin/env bash
# Rebuild and rerun the second fractional-design instance from the claimgate
# v0.1.0 tag.
#
# The script is read-only with respect to the claimgate repository: it only
# reads git metadata and pipes `git archive` out of it. Everything it writes
# lives under this record's directory.
#
# Usage:
#   bash reproduce.sh
#
# Environment:
#   CLAIMGATE_REPO  path to the claimgate clone (default: ../claimgate beside
#                   the repository that holds this record)
#   CLAIMGATE_TAG   tag to build against (default v0.1.0)
#
# Requirements: Git, CMake 3.20 or newer, a C++17 compiler, and the claimgate
# v0.1.0 tag in the clone.
#
# Every command and its exit code is appended to commands.txt, with paths
# written relative to this directory. Each command's output is kept in logs
# (the instance run's output is also results/run-all.log). The script exits
# non-zero if any command failed.

set -u -o pipefail

RECORD_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$RECORD_DIR"
HARNESS_REPO="$(cd "$RECORD_DIR/../.." && pwd)"
CLAIMGATE_REPO="${CLAIMGATE_REPO:-$(cd "$HARNESS_REPO/.." && pwd)/claimgate}"
CLAIMGATE_TAG="${CLAIMGATE_TAG:-v0.1.0}"

SNAP="$ROOT/src-claimgate-0.1.0"
BUILD="$ROOT/build"
INSTALL="$ROOT/install/claimgate-0.1.0"
RESULTS="$ROOT/results"
LOGS="$ROOT/logs"
COMMANDS="$ROOT/commands.txt"
EXE="$BUILD/instance/Release/instance.exe"

mkdir -p "$ROOT" "$LOGS"

rel() {
    printf '%s\n' "$1" | sed "s|$ROOT|.|g; s|$CLAIMGATE_REPO|\$CLAIMGATE_REPO|g"
}

FAILURES=0
COMMAND_COUNT=0

{
    echo "d-second-instance record: commands and exit codes"
    echo "date: $(date -u '+%Y-%m-%dT%H:%M:%SZ')"
    echo "host: $(uname -srm)"
    echo "claimgate repo: $(rel "$CLAIMGATE_REPO")"
    echo "claimgate tag: $CLAIMGATE_TAG"
    echo "instance root: $(rel "$ROOT")"
    echo
} > "$COMMANDS"

run() {
    local name="$1"
    local log="$2"
    local command="$3"

    COMMAND_COUNT=$((COMMAND_COUNT + 1))
    {
        echo "[$COMMAND_COUNT] $name"
        echo "  command: $(rel "$command")"
        echo "  cwd: $(rel "$ROOT")"
    } >> "$COMMANDS"

    ( cd "$ROOT" && bash -c "$command" ) > "$log" 2>&1
    local status=$?

    echo "  exit: $status" >> "$COMMANDS"
    echo >> "$COMMANDS"

    if [ "$status" -ne 0 ]; then
        FAILURES=$((FAILURES + 1))
        echo "FAILED ($status): $command" >&2
        echo "--- tail of $log ---" >&2
        tail -n 40 "$log" >&2
        echo "--- end of tail ---" >&2
    fi
}

# --- tool versions and the repository state ---------------------------------

run version-git "$LOGS/version-git.log" "git --version"
run version-cmake "$LOGS/version-cmake.log" "cmake --version"
run claimgate-status "$LOGS/claimgate-status.log" \
    "git -C '$CLAIMGATE_REPO' status --porcelain"
run claimgate-head "$LOGS/claimgate-head.log" \
    "git -C '$CLAIMGATE_REPO' rev-parse HEAD"
run claimgate-tag-commit "$LOGS/claimgate-tag-commit.log" \
    "git -C '$CLAIMGATE_REPO' rev-parse '$CLAIMGATE_TAG^{commit}'"
run claimgate-tag-object "$LOGS/claimgate-tag-object.log" \
    "git -C '$CLAIMGATE_REPO' rev-parse '$CLAIMGATE_TAG'"

# --- clean the generated directories ----------------------------------------

run clean-snapshot "$LOGS/clean-snapshot.log" "rm -rf '$SNAP'"
run clean-build "$LOGS/clean-build.log" "rm -rf '$BUILD' '$INSTALL' '$RESULTS'"
run make-dirs "$LOGS/make-dirs.log" "mkdir -p '$SNAP' '$BUILD' '$INSTALL' '$RESULTS'"
run extract-tag "$LOGS/extract-tag.log" \
    "git -C '$CLAIMGATE_REPO' archive '$CLAIMGATE_TAG' | tar -x -C '$SNAP'"

# --- build, test and install claimgate 0.1.0 from the tag -------------------

run claimgate-configure "$LOGS/claimgate-configure.log" \
    "cmake -S '$SNAP' -B '$BUILD/claimgate-msvc' -DCMAKE_BUILD_TYPE=Release"
run claimgate-build "$LOGS/claimgate-build.log" \
    "cmake --build '$BUILD/claimgate-msvc' --config Release"
run claimgate-ctest "$LOGS/claimgate-ctest.log" \
    "ctest --test-dir '$BUILD/claimgate-msvc' -C Release --output-on-failure"
run claimgate-example "$LOGS/claimgate-example.log" \
    "'$BUILD/claimgate-msvc/Release/claimgate-example.exe'"
run claimgate-install "$LOGS/claimgate-install.log" \
    "cmake --install '$BUILD/claimgate-msvc' --config Release --prefix '$INSTALL'"

# --- build the instance against the installed package -----------------------

run instance-configure "$LOGS/instance-configure.log" \
    "cmake -S '$RECORD_DIR' -B '$BUILD/instance' -DCMAKE_BUILD_TYPE=Release -DCMAKE_PREFIX_PATH='$INSTALL'"
run instance-build "$LOGS/instance-build.log" \
    "cmake --build '$BUILD/instance' --config Release"

# --- run the full factorial, the fraction, and the certificate --------------
# The run log is the record's primary result file, so it is written directly
# into results/ as well as being the command log.

run instance-run "$RESULTS/run-all.log" "'$EXE' all --out '$RESULTS'"

# --- confirm the claimgate repository was not modified ----------------------

run claimgate-status-after "$LOGS/claimgate-status-after.log" \
    "git -C '$CLAIMGATE_REPO' status --porcelain"

# --- list what the run produced ---------------------------------------------

run results-listing "$LOGS/results-listing.log" \
    "find '$RESULTS' -type f -printf '%P %s bytes\\n' | sort"

{
    echo "total commands: $COMMAND_COUNT"
    echo "failures: $FAILURES"
} >> "$COMMANDS"

echo "commands recorded in $COMMANDS ($COMMAND_COUNT commands, $FAILURES failures)"

if [ "$FAILURES" -ne 0 ]; then
    exit 1
fi
exit 0
