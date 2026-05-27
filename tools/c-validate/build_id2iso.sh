#!/usr/bin/env bash
# Build the id2iso cross-validation harness.
# Only requires mp.c and its dependencies (no quaternion, no signing).
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
REF="$REPO_ROOT/reference"
BUILD="$REF/build_sign"

if [ ! -f "$BUILD/src/gf/ref/lvl1/libsqisign_gf_lvl1.a" ]; then
    echo "ERROR: C reference not built. Run:" >&2
    echo "  cd $REF && mkdir -p build_sign && cd build_sign && cmake .. -DGMP_LIBRARY=MINI -DENABLE_SIGN=ON -DCMAKE_BUILD_TYPE=Release -DSQISIGN_BUILD_TYPE=ref && make -j\$(nproc)" >&2
    exit 1
fi

# Include paths
INC_COMMON="$REF/src/common/generic/include"
INC_NAMESPACE="$REF/include"
INC_PRECOMP="$REF/src/precomp/ref/lvl1/include"
INC_GF="$REF/src/gf/ref/include"
INC_MP="$REF/src/mp/ref/generic/include"

# Source: compile mp.c directly
MP_SRC="$REF/src/mp/ref/generic/mp.c"

# Libraries needed for mp.c (field arithmetic)
GF_LIB="$BUILD/src/gf/ref/lvl1/libsqisign_gf_lvl1.a"
PRECOMP_LIB="$BUILD/src/precomp/ref/lvl1/libsqisign_precomp_lvl1.a"

cc -O2 -o "$SCRIPT_DIR/id2iso_cval" \
    "$SCRIPT_DIR/id2iso_validate.c" \
    "$MP_SRC" \
    -DRADIX_64 \
    -DDISABLE_NAMESPACING \
    -I"$INC_COMMON" \
    -I"$INC_NAMESPACE" \
    -I"$INC_PRECOMP" \
    -I"$INC_GF" \
    -I"$INC_MP" \
    "$GF_LIB" \
    "$PRECOMP_LIB" \
    -lm

echo "Built: $SCRIPT_DIR/id2iso_cval"
