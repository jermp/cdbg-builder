#!/usr/bin/env bash
# End-to-end test: build the tiny dataset and verify cdgb-build's output.
# Run from anywhere; binary is expected at <repo>/build/cdgb-build.

set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
REPO="$(cd "$HERE/.." && pwd)"
BIN="$REPO/build/cdgb-build"
K=5

if [[ ! -x "$BIN" ]]; then
    echo "error: $BIN not found. Build the project first:" >&2
    echo "  mkdir -p $REPO/build && cd $REPO/build && cmake .. && make -j" >&2
    exit 2
fi

WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT

FILELIST="$WORK/filenames.txt"
printf '%s\n' \
    "$HERE/inputs/a.fa" \
    "$HERE/inputs/b.fa" \
    "$HERE/inputs/c.fa" \
    > "$FILELIST"

OUT="$WORK/tiny"
"$BIN" -i "$FILELIST" -k "$K" -o "$OUT" -t 1

python3 "$HERE/verify.py" --filenames "$FILELIST" --out "$OUT" -k "$K"
