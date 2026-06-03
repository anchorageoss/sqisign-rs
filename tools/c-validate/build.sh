#!/bin/bash
# Build the C cross-validation harness.
#
# Extracts the static helpers (lines 1..518) from
# reference/src/gf/ref/lvl1/fp_p5248_64.c so we can compile the
# self-contained portion without GMP / the reference build system,
# then compiles cval.c which includes that slice.

set -euo pipefail
cd "$(dirname "$0")"

REF_FP="../../reference/src/gf/ref/lvl1/fp_p5248_64.c"
if [ ! -f "$REF_FP" ]; then
    echo "error: $REF_FP not found; run \`git submodule update --init\` first" >&2
    exit 1
fi

# Stop at the line that pulls in fp.h (the public API wrappers start
# right after).
awk '/^#include <fp.h>/{exit} {print}' "$REF_FP" > fp_p5248_64_static.c

cc -O2 -Wno-unused-function -o cval cval.c
echo "built ./cval"
