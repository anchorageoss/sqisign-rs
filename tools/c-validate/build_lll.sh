#!/usr/bin/env bash
# Build the LLL/normeq/lat_ball cross-validation harness.
# Requires: the C reference built with ENABLE_SIGN=ON in reference/build_sign/
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
REF="$REPO_ROOT/reference"
BUILD="$REF/build_sign"

if [ ! -f "$BUILD/src/quaternion/ref/generic/libsqisign_quaternion_generic.a" ]; then
    echo "ERROR: C reference not built with signing. Run:" >&2
    echo "  cd $REF && mkdir -p build_sign && cd build_sign && cmake .. -DGMP_LIBRARY=MINI -DENABLE_SIGN=ON -DCMAKE_BUILD_TYPE=Release -DSQISIGN_BUILD_TYPE=ref && make -j\$(nproc)" >&2
    exit 1
fi

# Include paths
INC_COMMON="$REF/src/common/generic/include"
INC_INTBIG="$REF/src/quaternion/ref/generic/include"
INC_INTBIG_INT="$REF/src/quaternion/ref/generic/internal_quaternion_headers"
INC_MINIGMP="$REF/src/mini-gmp"
INC_NAMESPACE="$REF/include"

# Libraries and objects
QUAT_LIB="$BUILD/src/quaternion/ref/generic/libsqisign_quaternion_generic.a"
MINIGMP_OBJ="$BUILD/CMakeFiles/GMP.dir/src/mini-gmp/mini-gmp.c.o"
MINIGMP_EXTRA_OBJ="$BUILD/CMakeFiles/GMP.dir/src/mini-gmp/mini-gmp-extra.c.o"

# Find the RNG and crypto objects
RNG_OBJ="$BUILD/src/common/generic/CMakeFiles/sqisign_common_sys.dir/randombytes_system.c.o"
FIPS_OBJ="$BUILD/src/common/generic/CMakeFiles/sqisign_common_sys.dir/fips202.c.o"
AES_OBJ=""
if [ -f "$BUILD/src/common/generic/CMakeFiles/sqisign_common_test.dir/__/ref/aes_c.c.o" ]; then
    AES_OBJ="$BUILD/src/common/generic/CMakeFiles/sqisign_common_test.dir/__/ref/aes_c.c.o"
fi

cc -O2 -o "$SCRIPT_DIR/lll_cval" \
    "$SCRIPT_DIR/lll_validate.c" \
    -DMINI_GMP \
    -DRADIX_64 \
    -I"$INC_COMMON" \
    -I"$INC_INTBIG" \
    -I"$INC_INTBIG_INT" \
    -I"$INC_MINIGMP" \
    -I"$INC_NAMESPACE" \
    "$QUAT_LIB" \
    "$MINIGMP_OBJ" \
    "$MINIGMP_EXTRA_OBJ" \
    "$RNG_OBJ" \
    "$FIPS_OBJ" \
    $AES_OBJ \
    -lm

echo "Built: $SCRIPT_DIR/lll_cval"
