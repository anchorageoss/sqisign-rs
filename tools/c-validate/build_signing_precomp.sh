#!/usr/bin/env bash
# Build the signing precomputed constants cross-validation harness.
# Requires: the C reference built with ENABLE_SIGN=ON in reference/build_sign/
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
REF="$REPO_ROOT/reference"
BUILD="$REF/build_sign"

if [ ! -f "$BUILD/src/precomp/ref/lvl1/libsqisign_precomp_lvl1.a" ]; then
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
INC_PRECOMP="$REF/src/precomp/ref/lvl1/include"
INC_GF="$REF/src/gf/ref/include"
INC_EC="$REF/src/ec/ref/include"
INC_HD="$REF/src/hd/ref/include"

# Libraries
PRECOMP_LIB="$BUILD/src/precomp/ref/lvl1/libsqisign_precomp_lvl1.a"
QUAT_LIB="$BUILD/src/quaternion/ref/generic/libsqisign_quaternion_generic.a"
GF_LIB="$BUILD/src/gf/ref/lvl1/libsqisign_gf_lvl1.a"
EC_LIB="$BUILD/src/ec/ref/lvl1/libsqisign_ec_lvl1.a"
MINIGMP_OBJ="$BUILD/CMakeFiles/GMP.dir/src/mini-gmp/mini-gmp.c.o"
MINIGMP_EXTRA_OBJ="$BUILD/CMakeFiles/GMP.dir/src/mini-gmp/mini-gmp-extra.c.o"

# RNG and crypto objects
RNG_OBJ="$BUILD/src/common/generic/CMakeFiles/sqisign_common_sys.dir/randombytes_system.c.o"
FIPS_OBJ="$BUILD/src/common/generic/CMakeFiles/sqisign_common_sys.dir/fips202.c.o"

cc -O2 -o "$SCRIPT_DIR/signing_precomp_cval" \
    "$SCRIPT_DIR/signing_precomp_validate.c" \
    -DMINI_GMP \
    -DRADIX_64 \
    -DDISABLE_NAMESPACING \
    -I"$INC_COMMON" \
    -I"$INC_INTBIG" \
    -I"$INC_INTBIG_INT" \
    -I"$INC_MINIGMP" \
    -I"$INC_NAMESPACE" \
    -I"$INC_PRECOMP" \
    -I"$INC_GF" \
    -I"$INC_EC" \
    -I"$INC_HD" \
    "$PRECOMP_LIB" \
    "$GF_LIB" \
    "$EC_LIB" \
    "$QUAT_LIB" \
    "$MINIGMP_OBJ" \
    "$MINIGMP_EXTRA_OBJ" \
    "$RNG_OBJ" \
    "$FIPS_OBJ" \
    -lm

echo "Built: $SCRIPT_DIR/signing_precomp_cval"
