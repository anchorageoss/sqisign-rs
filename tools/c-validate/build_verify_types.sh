#!/bin/bash
# Build the verification types cross-validation harness.

set -euo pipefail
cd "$(dirname "$0")"

REF="../../reference"

if [ ! -f "$REF/src/verification/ref/lvlx/common.c" ]; then
    echo "error: reference submodule not found; run \`git submodule update --init\` first" >&2
    exit 1
fi

# Copy source files
cp -f "$REF/src/gf/ref/lvl1/fp_p5248_64.c"  fp_p5248_64.c
cp -f "$REF/src/gf/ref/lvlx/fp.c"            fp_select.c
cp -f "$REF/src/gf/ref/lvlx/fp2.c"           fp2.c
cp -f "$REF/src/mp/ref/generic/mp.c"          mp.c
cp -f "$REF/src/ec/ref/lvlx/ec.c"            ec.c
cp -f "$REF/src/precomp/ref/lvl1/ec_params.c" ec_params.c
cp -f "$REF/src/common/generic/fips202.c"     fips202.c
cp -f "$REF/src/verification/ref/lvlx/encode_verification.c" encode_verification.c
cp -f "$REF/src/verification/ref/lvlx/common.c" common.c

cc -O2 -Wno-unused-function -Wno-unused-variable -Wno-unused-but-set-variable \
   -DDISABLE_NAMESPACING -DRADIX_64 -DHAVE_UINT128 -DNDEBUG \
   -I "$REF/include" \
   -I "$REF/src/common/generic/include" \
   -I "$REF/src/common/ref/include" \
   -I "$REF/src/precomp/ref/lvl1/include" \
   -I "$REF/src/gf/ref/include" \
   -I "$REF/src/mp/ref/generic/include" \
   -I "$REF/src/ec/ref/include" \
   -I "$REF/src/verification/ref/include" \
   -o verify_types_cval verify_types_validate.c

echo "built ./verify_types_cval"
