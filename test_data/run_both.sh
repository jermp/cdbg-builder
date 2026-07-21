#!/usr/bin/env bash
# Correctness + timing harness for a single dataset directory.
# Usage: run_both.sh <dataset_dir> <glob> [k] [threads] [extra cdbg-build args...]
# Prints verify.py result and keeps the build's verbose output for timing/RSS.
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
REPO="$(cd "$HERE/.." && pwd)"
BIN="$REPO/build/cdbg-build"
DIR="$1"; GLOB="$2"; K="${3:-31}"; THREADS="${4:-4}"
shift 4 || shift $#
EXTRA=("$@")

WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT
FILELIST="$WORK/filenames.txt"
# shellcheck disable=SC2086
ls $DIR/$GLOB | grep -v '/\._' > "$FILELIST"
echo "== dataset $DIR : $(wc -l < "$FILELIST") files, k=$K, t=$THREADS, extra=${EXTRA[*]:-none} =="
OUT="$WORK/out"
"$BIN" -i "$FILELIST" -k "$K" -o "$OUT" -t "$THREADS" -d "$WORK/buckets" --verbose "${EXTRA[@]}" 2>&1 \
  | grep -E '\[bucket-process\]|peak RSS|walk |load |resolve |total construction|num_color_sets' || true
python3 "$HERE/verify.py" --filenames "$FILELIST" --out "$OUT" -k "$K"
