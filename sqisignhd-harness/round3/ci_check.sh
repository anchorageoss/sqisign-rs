#!/usr/bin/env bash
# The compact format's one CI line: every self-generated level I vector is
# accepted by the SQIsignHD library's Sage verifier at our parameters, and
# every tampered one rejected. Pass/fail only. Needs `sage` on the path (or
# SAGE=/path/to/sage) and network access to clone the library once.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SAGE="${SAGE:-sage}"
LIB="${SQISIGNHD_LIB:-$HERE/SQISignHD-lib}"
if [ ! -d "$LIB/Verification" ]; then
  git clone -q https://github.com/Pierrick-Dartois/SQISignHD-lib "$LIB"
  git -C "$LIB" config submodule.Verification/Theta_dim4.url https://github.com/Pierrick-Dartois/Theta_dim4
  git -C "$LIB" submodule update -q --init
fi
"$SAGE" -python "$HERE/verify_r3.py" --lib "$LIB" --pk "$HERE/Data/Public_keys_r3lvl1.txt" --sig "$HERE/Data/Signatures_r3lvl1.txt" -n 5
"$SAGE" -python "$HERE/verify_r3.py" --lib "$LIB" --pk "$HERE/Data/Public_keys_r3lvl1.txt" --sig "$HERE/Data/Signatures_r3lvl1_bad.txt" -n 5 --expect-reject
echo "compact vectors: accepted set accepted, tampered set rejected"
