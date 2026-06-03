#!/bin/bash
# Build the basis cross-validation harness.

set -euo pipefail
cd "$(dirname "$0")"

REF="../../reference"

if [ ! -f "$REF/src/ec/ref/lvlx/ec.c" ]; then
    echo "error: reference submodule not found; run \`git submodule update --init\` first" >&2
    exit 1
fi

# Copy source files
cp "$REF/src/gf/ref/lvl1/fp_p5248_64.c"  fp_p5248_64.c
cp "$REF/src/gf/ref/lvlx/fp.c"            fp_select.c
cp "$REF/src/gf/ref/lvlx/fp2.c"           fp2.c
cp "$REF/src/mp/ref/generic/mp.c"          mp.c
cp "$REF/src/ec/ref/lvlx/ec.c"            ec.c
cp "$REF/src/ec/ref/lvlx/ec_jac.c"        ec_jac.c
cp "$REF/src/ec/ref/lvlx/xisog.c"         xisog.c
cp "$REF/src/ec/ref/lvlx/xeval.c"         xeval.c
cp "$REF/src/ec/ref/lvlx/isog_chains.c"   isog_chains.c
cp "$REF/src/ec/ref/lvlx/basis.c"          basis.c
cp "$REF/src/precomp/ref/lvl1/ec_params.c" ec_params.c
cp "$REF/src/precomp/ref/lvl1/e0_basis.c"  e0_basis.c

cc -O2 -Wno-unused-function -Wno-unused-variable -Wno-unused-but-set-variable \
   -DDISABLE_NAMESPACING -DRADIX_64 -DHAVE_UINT128 -DNDEBUG \
   -I "$REF/include" \
   -I "$REF/src/common/generic/include" \
   -I "$REF/src/precomp/ref/lvl1/include" \
   -I "$REF/src/gf/ref/include" \
   -I "$REF/src/mp/ref/generic/include" \
   -I "$REF/src/ec/ref/include" \
   -o basis_cval basis_validate.c

echo "built ./basis_cval"
