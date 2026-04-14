#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DIAG_DIR="$ROOT_DIR/diagrams"

for dot in "$DIAG_DIR"/*.dot; do
  base="${dot%.dot}"
  dot -Tpng "$dot" -o "${base}.png"
  dot -Tsvg "$dot" -o "${base}.svg"
  echo "Generated: ${base}.png and ${base}.svg"
done

uv run python "$DIAG_DIR/platform_visual_diagrams.py"
echo "Generated: platform_visual_diagrams.png"
