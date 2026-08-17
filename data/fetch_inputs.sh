#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# fetch_inputs.sh — retrieve the AdK alignments and reference structures.
#
# These files come from Guan et al. (2024), AFF2024. That repository carries no
# licence file, so the inputs are fetched here rather than redistributed.
# Checksums are verified so you can confirm you have the exact files this
# analysis used.
#
# Usage:  bash data/fetch_inputs.sh
# Result: data/AFF2024/AdK/...
# ---------------------------------------------------------------------------
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DEST="$HERE/AFF2024"
UPSTREAM="https://github.com/Gxy-with-luv/AFF2024"

if [ -d "$DEST/.git" ]; then
  echo "==> $DEST already exists, skipping clone"
else
  echo "==> cloning $UPSTREAM (shallow)"
  git clone --depth 1 "$UPSTREAM" "$DEST"
fi

echo "==> verifying checksums"
cd "$DEST/AdK"
if command -v shasum >/dev/null 2>&1; then
  shasum -a 256 -c "$HERE/input_checksums.sha256"
else
  sha256sum -c "$HERE/input_checksums.sha256"
fi

echo
echo "==> inputs ready:"
echo "    $DEST/AdK/Frustration-Energy-space/MSA.a3m        (2651 x 214)"
echo "    $DEST/AdK/AdK-mixing/Closed500.a3m                (501 seqs)"
echo "    $DEST/AdK/AdK-mixing/Open100.a3m                  (101 seqs)"
echo "    $DEST/AdK/Frustration-Energy-space/{closed,open}CA.pdb"
