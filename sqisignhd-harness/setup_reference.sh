#!/usr/bin/env bash
# setup_reference.sh -- fetch the SQIsignHD reference verifier + Theta_dim4
# submodule into a local checkout that extract_vectors.py can use.
#
# The upstream .gitmodules pins Theta_dim4 to an *ssh* URL
# (git@github.com:Pierrick-Dartois/Theta_dim4.git), which fails for anonymous
# clones, so we clone the submodule explicitly over https instead of using
# `git submodule update`.
#
# Usage:  ./setup_reference.sh [DEST]
#   DEST defaults to ./SQISignHD-lib  (or $SQISIGNHD_LIB if set).
set -euo pipefail

DEST="${1:-${SQISIGNHD_LIB:-$(pwd)/SQISignHD-lib}}"
MAIN_URL="https://github.com/Pierrick-Dartois/SQISignHD-lib.git"
THETA_URL="https://github.com/Pierrick-Dartois/Theta_dim4.git"

echo "Destination: $DEST"

if [ ! -d "$DEST/.git" ]; then
  echo "Cloning $MAIN_URL ..."
  git clone --depth 1 "$MAIN_URL" "$DEST"
else
  echo "Main repo already present, skipping clone."
fi

THETA_DIR="$DEST/Verification/Theta_dim4"
if [ ! -d "$THETA_DIR/Theta_dim4_sage" ]; then
  echo "Cloning Theta_dim4 submodule (https) ..."
  rm -rf "$THETA_DIR"
  git clone --depth 1 "$THETA_URL" "$THETA_DIR"
else
  echo "Theta_dim4 submodule already present, skipping clone."
fi

echo
echo "Reference ready at: $DEST"
echo "Main commit:    $(git -C "$DEST" rev-parse HEAD)"
echo "Theta_dim4:     $(git -C "$THETA_DIR" rev-parse HEAD)"
echo
echo "Next:"
echo "  export SQISIGNHD_LIB=\"$DEST\""
echo "  sage extract_vectors.py --lvl 1 --n 5 --out test_vectors_l1.json"
echo "  python3 validate_vectors.py test_vectors_l1.json"
