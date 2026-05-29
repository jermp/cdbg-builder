#!/usr/bin/env bash
# Run all per-phase test executables. Expects them built in <repo>/build
# (cmake .. && make -j). Returns non-zero if any test fails.

set -uo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
REPO="$(cd "$HERE/.." && pwd)"
BUILD="$REPO/build"

tests=(
    test_stitch
)

rc=0
for t in "${tests[@]}"; do
    bin="$BUILD/$t"
    if [[ ! -x "$bin" ]]; then
        echo "error: $bin not found; build first (cd build && cmake .. && make -j)" >&2
        rc=2
        continue
    fi
    echo "=== $t ==="
    if ! "$bin"; then
        echo "  -> $t FAILED" >&2
        rc=1
    fi
done

if [[ $rc -eq 0 ]]; then
    echo "all tests passed"
fi
exit $rc
