# All-zero vectors in the dimension-4 verifier, 2026-11-16

Counts taken on the code of tag `paper-dim4-2026-11` with the counters of
`crates/verify/src/hd/zero_trace.rs` added (feature `instrument`, no change
of behaviour), before the fixes; harness `paper-dim4/bench/src/zero_trace.rs`.

- `at-tag/level{1,3,5}.txt`: 20, 8 and 4 signatures. Every bit of the
  signature, the key and the message; every bit of either curve with the
  canonical hints of the new curve; the commitment curve with the hint check
  off; top bits and negations of the scalars; `a` and `b` both even; scalars
  of zeros and ones; random scalars; other good degrees, among them
  `2^e - p` for small `p` and for `v2(a2)` up to and past the bound 40;
  crafted curves (`A` = 0, 6, 2, -2, i, the other curve, conjugates). Each
  input in the three verifier modes.
- `before-after/`: a second set of inputs (3, 2 and 1 signatures) run before
  and after the fixes. The records are identical.
- `unit-tests-at-tag.txt`: the helpers called directly with all-zero inputs.
- `timing/`: two pinned passes of the fixed code against the clean session.
- Tests on the fixed code: 161 passed, 0 failed, 3 ignored over the workspace
  at the three levels (known-answer, tampered and compact conformance tests
  among them). The SQIsignHD-lib cross-check at level I was not run: the
  library is no longer on this machine and fetching it was not permitted.

## Result

159 258 tampered inputs, 476 776 runs: no all-zero theta vector, null point,
Kummer point or Jacobian triple at any of the 28 sites, in any run. 222 runs
accepted, all of them the two malleabilities the paper states: the full
negation (32 inputs, every one accepted in every mode) and the top bit of
`a` or `b`, alone or with other top bits (the acceptances the harness marks
UNEXPECTED in `sig-bit` are those same two bits reached through the bytes:
bits 912 and 1001, 1405 and 1539, 2031). The top bit of `c_or_d` alone was
never accepted. Nothing else was accepted.

## Sites fixed

| Site | Before | After |
|---|---|---|
| `point.rs` `projective_eq` (callers: `middle_codomain_matches`, `gluing_images_paired`) | true when either side is all zero | false |
| `point.rs` `normalize` | inverts a zero pivot | `None` |
| `two_point.rs` `null_fixed_by_pi` | true on the zero vector | false |
| `self_contained.rs` `x_eq` | true for `(0, y, 0)` against anything, and for two identities | false when either `z = 0` |
| `self_contained.rs` image check, fourth component | `z = 0` | `z = 0`, `x` and `y` not zero |
| `self_contained.rs` `dim1_null` | inverts `z` untested | `None` on `z = 0`, a reject |
| `dim2.rs` `dual_eval`, last case | the pair of identities for the zero vector and for non-product coordinates with `Z1Z2 = 0` | `None` |
| `dim2.rs` `lift_kummer` | `(0 : 0)` read as the identity | `None` |
| `dim2.rs` `tuple_product_xz` | a Jacobian `(0, y, 0)` read as the identity | `None` |
| `dim2.rs` `FoldedBaseChangeDim2` inverse | products of the others of `zero12` with a zero among them | `None`, the dense path decides |
| `gluing.rs` `special_image` | could return the zero vector | `None` |
| `product_theta.rs` `montgomery_to_theta` (no caller) | `(0 : 0)` read as the identity | `None` |

## Sites reported, not changed

- `basis.rs` `jac_to_affine`: inverts `z^2` and `z^3` without a test. Its
  one caller on the verification path now tests `z`.
- `response.rs` `recover_response_prof`: inverts `w_chal`. The verifier
  checks that `w_chal` has exact order `2^r` first; a direct caller of the
  public function gets no test.
- `dual_eval` with `Z1Z2` not zero, and `product_to_theta_points_dim4_dim2`:
  neither tests that its input is a product point. Not a zero pattern.
- `ec/jacobian.rs` (`to_xz`, `jac_add`, `jac_dbl`): the triple `(0, y, 0)`,
  the all-zero one included, is read as the identity. This is the round-3
  dimension-2 code under the known-answer tests; the dimension-4 layer now
  rejects such a triple where it enters.
- Every other inversion in `hd/` follows a zero test or has nonzero input by
  construction: `field::inv` in `dim4.rs`, `mat_inverse`,
  `special_compute_codomain`, `dual_eval`; `batched_inv` in `isogeny.rs` and
  `dim2.rs`; the products of the others in `structure.rs`, `eval_f2_dual`,
  `splitting_eval`, `dual_eval`; `inv_odd` on `a` or `b` (not both even), on
  `q` (odd for a good degree), on `a1` (odd) and on unit pivots in `kani.rs`.
- `double_assign` and `apply_plain_image` send the zero vector to itself
  without a signal; their inputs are nonzero now that the gluing image is.
