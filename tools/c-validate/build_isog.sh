#!/bin/bash
# Build the isogeny cross-validation harness.
#
# Copies the needed C source files from the reference and compiles
# isog_validate.c which #includes them all in a single translation unit.

set -euo pipefail
cd "$(dirname "$0")"

REF="../../reference"

# Verify submodule is present
if [ ! -f "$REF/src/ec/ref/lvlx/ec.c" ]; then
    echo "error: reference submodule not found; run \`git submodule update --init\` first" >&2
    exit 1
fi

# Copy source files we need
cp "$REF/src/gf/ref/lvl1/fp_p5248_64.c"  fp_p5248_64.c
cp "$REF/src/gf/ref/lvlx/fp.c"            fp_select.c
cp "$REF/src/gf/ref/lvlx/fp2.c"           fp2.c
cp "$REF/src/mp/ref/generic/mp.c"          mp.c
cp "$REF/src/ec/ref/lvlx/ec.c"            ec.c
cp "$REF/src/ec/ref/lvlx/ec_jac.c"        ec_jac.c
cp "$REF/src/ec/ref/lvlx/xisog.c"         xisog.c
cp "$REF/src/ec/ref/lvlx/xeval.c"         xeval.c
cp "$REF/src/ec/ref/lvlx/isog_chains.c"   isog_chains.c

cc -O2 -Wno-unused-function -Wno-unused-variable -Wno-unused-but-set-variable \
   -DDISABLE_NAMESPACING -DRADIX_64 -DHAVE_UINT128 -DNDEBUG \
   -I "$REF/include" \
   -I "$REF/src/common/generic/include" \
   -I "$REF/src/precomp/ref/lvl1/include" \
   -I "$REF/src/gf/ref/include" \
   -I "$REF/src/mp/ref/generic/include" \
   -I "$REF/src/ec/ref/include" \
   -o isog_cval isog_validate.c

echo "built ./isog_cval"
