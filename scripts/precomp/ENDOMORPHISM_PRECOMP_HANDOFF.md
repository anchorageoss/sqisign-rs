# Endomorphism Precomputation: Status and Handoff

## Goal

Generate `endomorphism_action.rs` for all 3 security levels (1, 3, 5) using pure SageMath
instead of parsing C reference output. This eliminates the C build dependency and ensures
the Rust constants are generated directly from the mathematical specification.

## Current Status: BLOCKED ON COMPUTE TIME

The code is ready and the first q-value (q=1) succeeds for level 1. The q=5 computation
also succeeds through the qlapoti phase (~15 min) but the full pipeline for all 7 q-values
per level will take 1-3+ hours. This needs to run on a GCP VM overnight.

## What Works

1. **e0_basis generation** (COMPLETE): `generate_e0_basis_rust.sage` generates e0_basis.rs
   for all 3 levels, byte-identical to C reference output.

2. **qlapoti fix** (COMPLETE): Fixed Singular overflow bug in
   `scripts/precomp/deps/deuring-2D/qlapoti/qlapoti.py` (lines 48-62).
   - Original code: `Zmod(N)` polynomial ring triggers Singular, which uses C long internally
     and overflows for N > 2^63 (~250-bit SQIsign primes)
   - Fix: Extract polynomial coefficients mod N via numerical evaluation at 4 points,
     using `qq_mod_N()` helper that handles QQ->ZZ conversion via modular inversion
   - Self-test passes (exit 0) on random ~100-bit primes
   - Real usage: qlapoti completes successfully for level 1 q=5 (~15 min compute)

3. **deuring2d field extension fix** (COMPLETE): Fixed TypeError in
   `run_endomorphism_precomp.sage` where `ctx.iota` and `ctx.pi` were defined over Fp2
   but `ctx.E0` was changed to Fbig (extension field).
   - Fix: Re-derive `ctx.iota = E0_big.automorphisms()[-1]` and
     `ctx.pi = E0_big.frobenius_isogeny()` after changing `ctx.E0`

4. **Runner script**: `run_all_endomorphism.sh` runs all 3 levels sequentially with logging.

## What Needs To Happen

### Step 1: Run the computation (GCP VM, overnight)

```bash
cd /path/to/sqisign-rs

# If using micromamba (as on this machine):
export PATH="/path/to/sqisign-rs/bin:$PATH"
eval "$(micromamba shell hook --shell bash)"
micromamba activate sage

# Run all levels
bash scripts/precomp/run_all_endomorphism.sh
```

Expected output: `crates/precomp/src/level{1,3,5}/endomorphism_action.rs`

### Step 2: After computation succeeds

1. **Verify the generated files compile**:
   ```bash
   cargo build -p sqisign-precomp --features signing
   ```

2. **Run the existing validation tests**:
   ```bash
   cargo test -p sqisign-precomp --features signing
   ```

3. **Compare against C-parser-generated constants** (the current fallback):
   The C-parser-generated endomorphism_action.rs files are in git history. Compare the
   Fp limb arrays to verify byte-identical output. The action matrices (BigInt values)
   should also match.

4. **If all tests pass**: commit the Sage-generated files and remove the C-parser dependency.

## Architecture

### File Layout

```
scripts/precomp/
  generate_e0_basis_rust.sage       # Working: generates e0_basis.rs
  run_endomorphism_precomp.sage     # Main script: generates endomorphism_action.rs
  run_all_endomorphism.sh           # Runner for all 3 levels
  parse_fp_constants.py             # FALLBACK: C parser (still works, used for current constants)
  deps/
    deuring-2D/                     # Cloned: github.com/ThetaIsogenies/deuring-2D
      qlapoti/
        qlapoti.py                  # PATCHED: Singular overflow fix
    two-isogenies/                  # Cloned: github.com/ThetaIsogenies/two-isogenies
      Theta-SageMath/               # Needed by deuring2d
```

### How the computation works

For each security level, there are N curves (7 for level 1) with known endomorphism rings.
Each curve is identified by a small prime q (1, 5, 17, 37, 41, 53, 97 for level 1).

For each curve:
1. **q=1**: Shortcut via automorphisms of E0 (fast, ~seconds)
2. **q>1**: Full Deuring correspondence computation:
   a. Construct quaternion ideal I from `maxorders.py` data
   b. Create Deuring2D context, possibly over field extension Fbig
   c. Call `ctx.IdealToIsogeny(I)` which internally:
      - Calls `qlapoti.solve()` to decompose the ideal (~15 min per q)
      - Computes 2-dimensional isogeny via theta functions
   d. Convert result to Montgomery form
   e. Compute endomorphism action matrices

### Montgomery form encoding

Critical detail: SQIsign uses non-standard radix for limb encoding:
- Level 1: 51-bit limbs, 5 words, R = 2^(51*5) = 2^255
- Level 3: 55-bit limbs, 7 words, R = 2^(55*7) = 2^385
- Level 5: 57-bit limbs, 9 words, R = 2^(57*9) = 2^513

This is NOT the standard 64-bit limb Montgomery form. The radix maps come from the C
reference's `cformat.py`.

### Output format

Each `endomorphism_action.rs` contains per-curve:
- 8 Fp limb arrays: curve coefficients (A, C=1, A24=(A+2)/4)
- 12 Fp limb arrays: torsion basis (P, Q, P-Q in projective X:Z form)
- 6 action matrices: [Lazy<BigInt>; 4] for i, j, k, gen2, gen3, gen4

## Known Issues / Risks

1. **Compute time**: Each q-value takes 5-30 minutes in qlapoti. With 7 q-values per level
   and 3 levels, total time is 2-10+ hours. Run overnight.

2. **iota selection**: The fix uses `E0_big.automorphisms()[-1]` to select the distortion
   map over the extended field. This mirrors the Deuring2D constructor's approach. If Sage
   orders automorphisms differently over Fbig vs Fp2, this could select the wrong one
   (the conjugate). If results are wrong, try `automorphisms()[2]` or filter by
   `scaling_factor()`.

3. **BinaryQF factoring**: The qlapoti inner loop calls `BinaryQF.solve_integer(rhs)` which
   needs to factor a ~250-bit number. This is the main bottleneck. If it gets stuck on a
   particular alpha, the code tries the next one (there's an outer loop over candidate alphas).

4. **Sage version**: Needs Sage >= 10.0. The `maxorders.py` from the SQIsign reference may
   need Sage >= 10.5 for some features. Current setup uses micromamba sage env.

5. **PARI memory**: Script allocates 16GB via `pari.allocatemem(1 << 34)`. The GCP VM needs
   at least 16GB RAM.

## Dependencies Setup (if cloning to a new machine)

```bash
cd scripts/precomp
mkdir -p deps
cd deps
git clone https://github.com/ThetaIsogenies/deuring-2D.git
cd deuring-2D && git submodule update --init && cd ..
git clone https://github.com/ThetaIsogenies/two-isogenies.git
cd ../..

# Apply the qlapoti fix (already applied in this repo):
# The patched file is scripts/precomp/deps/deuring-2D/qlapoti/qlapoti.py
```

## The qlapoti fix in detail

**File**: `scripts/precomp/deps/deuring-2D/qlapoti/qlapoti.py`, lines 48-75

**Original code** (crashes for N > 2^63):
```python
A,B = polygens(Zmod(N), 'A,B')
eqn1 = eqn.change_ring(Zmod(N))(A,B,0,0)
```

**Fixed code** (works for arbitrary N):
```python
def qq_mod_N(v):
    v = QQ(v)
    return ZZ(v.numerator()) * ZZ(v.denominator()).inverse_mod(N) % N
c00 = eqn(0,0,0,0)
c10 = eqn(1,0,0,0)
c01 = eqn(0,1,0,0)
c11 = eqn(1,1,0,0)
coeff_const = qq_mod_N(c00)
coeff_A = qq_mod_N(c10 - c00)
coeff_B = qq_mod_N(c01 - c00)
assert qq_mod_N(c11 - c10 - c01 + c00) == 0
```

**Why**: `Zmod(N)` for large N delegates to Singular's polynomial ring implementation,
which uses C `long` internally and overflows for N > 2^63. The fix avoids creating any
`Zmod(N)` polynomial ring by extracting the three needed coefficients (constant, coeff of A,
coeff of B) via numerical evaluation of the ZZ polynomial at (0,0,0,0), (1,0,0,0),
(0,1,0,0), and (1,1,0,0), then reducing mod N. The `qq_mod_N` helper handles the case where
coefficients are in QQ (from the `have1i2` substitution that introduces /2 denominators)
by computing `numerator * denominator^{-1} mod N`.

## Fallback

If the Sage computation fails or produces wrong results, the C-parser approach still works:

```bash
python3 scripts/precomp/parse_fp_constants.py
```

This parses the C reference's precomputed .c files directly. The current
`endomorphism_action.rs` files in git were generated this way and are known-good.
