# C cross-validation harnesses

Minimal C programs that compile parts of the C reference and run a
fixed sequence of operations on known inputs, printing each result as
hex-encoded bytes.  The Rust tests in each crate's `tests/` directory
perform the identical computations and assert byte-for-byte equality.

No GMP, no CMake, no reference build system required, each harness
compiles standalone from extracted source slices.

---

## sqisign-fp (field arithmetic)

**C harness:** `cval.c`
**Rust test:** `crates/fp/tests/c_crossvalidate.rs`

```
tools/c-validate/build.sh
tools/c-validate/cval
```

`build.sh` extracts the static helpers (lines 1..518) of
`reference/src/gf/ref/lvl1/fp_p5248_64.c` into
`fp_p5248_64_static.c`, then compiles `cval.c` against the slice.

Covers: `add`, `sub`, `neg`, `mul`, `sqr`, `inv`, `sqrt`.

---

## sqisign-ec (elliptic curve arithmetic)

**C harness:** `ec_validate.c`
**Rust test:** `crates/ec/tests/c_crossvalidate.rs`

```
tools/c-validate/build_ec.sh
tools/c-validate/ec_cval
```

`build_ec.sh` copies the needed C source files (fp, fp2, mp, ec,
ec\_jac) from the reference and compiles `ec_validate.c` as a single
translation unit that `#include`s them all.

Covers: `xDBL`, `xADD`, `xDBLADD`, `xMUL`, `ec_ladder3pt`,
`ec_j_inv` (two curves), `ec_normalize_point`, `jac_dbl`, `jac_add`,
`jac_to_xz`, `jac_add` with identity, `jac_add(P, -P)`.

---

## sqisign-ec isogeny layer

**C harness:** `isog_validate.c`
**Rust test:** `crates/ec/tests/c_crossvalidate_isog.rs`

```
tools/c-validate/build_isog.sh
tools/c-validate/isog_cval
```

`build_isog.sh` copies the needed C source files (fp, fp2, mp, ec,
ec\_jac, xisog, xeval, isog\_chains) from the reference and compiles
`isog_validate.c` as a single translation unit.

Covers: `xisog_2`, `xeval_2`, `xisog_4`, `xeval_4`, manual two-step
degree-2 chain, `ec_isomorphism`, `ec_iso_eval`, `ec_eval_small_chain`,
codomain A:C recovery, codomain j-invariant.

---

## sqisign-precomp (precomputed constants)

**C harness:** `precomp_validate.c`
**Rust test:** `crates/precomp/tests/c_crossvalidate_precomp.rs`

```
tools/c-validate/build_precomp.sh
tools/c-validate/precomp_cval
```

`build_precomp.sh` copies fp, fp2, mp, ec layers plus the Level 1
precomp files (ec\_params.c, e0\_basis.c) and compiles
`precomp_validate.c`.

Covers: scalar constants (TORSION\_EVEN\_POWER, cofactor),
E0 j-invariant, BASIS\_E0\_PX / BASIS\_E0\_QX encode/decode
round-trip, on-curve verification, cofactor clearing to 2^f-torsion.

---

## sqisign-ec basis generation

**C harness:** `basis_validate.c`
**Rust test:** `crates/ec/tests/c_crossvalidate_basis.rs`

```
tools/c-validate/build_basis.sh
tools/c-validate/basis_cval
```

`build_basis.sh` copies fp, fp2, mp, ec, ec\_jac, xisog, xeval,
isog\_chains, basis, ec\_params, and e0\_basis from the reference and
compiles `basis_validate.c`.

Covers: `ec_basis_E0_2f` (full and partial order), `ec_recover_y`
(two curves), `ec_curve_to_basis_2f_to_hint`, `ec_curve_to_basis_2f_from_hint`,
`is_on_curve`, `lift_basis`.

---

## sqisign-ec biextension (pairings + dlog)

**C harness:** `biext_validate.c`
**Rust test:** `crates/ec/tests/c_crossvalidate_biext.rs`

```
tools/c-validate/build_biext.sh
tools/c-validate/biext_cval
```

`build_biext.sh` copies the same sources as the basis harness plus
`biextension.c` and compiles `biext_validate.c`.

Covers: `weil`, `reduced_tate`, `clear_cofac`, `fp2_frob`,
`ec_dlog_2_weil`, `ec_dlog_2_tate` (full and partial torsion),
pairing order and bilinearity checks, dlog round-trip verification.

---

## Regenerating expected bytes

If the C reference changes (submodule update) and the Rust port
follows, rebuild and rerun the harness:

```
tools/c-validate/build_ec.sh && tools/c-validate/ec_cval
```

Paste the new hex strings into the corresponding Rust test file.  If
the Rust port drifts, the test fails with a hex diff identifying the
operation that disagrees.
