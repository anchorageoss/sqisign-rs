/*
 * Wrapper to compile the C reference biextension-test.c standalone.
 *
 * Build:  tools/c-validate/build_biext_reftest.sh
 * Run:    tools/c-validate/biext_reftest --seed=0
 */

#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <stdbool.h>
#include <time.h>
#include <assert.h>

/* Stub out tools.h timing functions */
clock_t tic(void) { return 0; }
float tac(void) { return 0; }
float TAC(const char *str) { (void)str; return 0; }
float toc(const clock_t t) { (void)t; return 0; }
float TOC(const clock_t t, const char *str) { (void)t; (void)str; return 0; }
float TOC_clock(const clock_t t, const char *str) { (void)t; (void)str; return 0; }
clock_t dclock(const clock_t t) { (void)t; return 0; }
float clock_to_time(const clock_t t, const char *str) { (void)t; (void)str; return 0; }
float clock_print(const clock_t t, const char *str) { (void)t; (void)str; return 0; }

/* Stub debug_print */
#define debug_print(msg) ((void)0)

/* Pull in the full fp layer */
#include "fp_p5248_64.c"
#include "fp_select.c"
#include "fp2.c"
#include "mp.c"

/* Provide the precomputed constants */
#include "ec_params.c"

/* Pull in EC + isogeny layers */
#include "ec.c"
#include "ec_jac.c"
#include "xisog.c"
#include "xeval.c"
#include "isog_chains.c"

/* Pull in basis generation */
#include "e0_basis.c"
#include "basis.c"

/* Pull in biextension (pairings + dlog) */
#include "biextension.c"

/* AES + CTR-DRBG RNG */
#include "aes_c.c"
#include "randombytes_ctrdrbg.c"

/* randombytes_select for seed generation (uses /dev/urandom) */
#include "randombytes_system.c"

/* Test extras */
#include "test_extras.c"

/* The actual test */
#include "../../reference/src/ec/ref/lvlx/test/biextension-test.c"
