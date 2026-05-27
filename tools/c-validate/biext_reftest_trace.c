/*
 * Instrumented version of the C reference biextension test.
 * Prints all intermediate values for cross-validation against Rust.
 *
 * Build: tools/c-validate/build_biext_reftest.sh  (but with this file)
 * Run:   tools/c-validate/biext_reftest_trace --seed=0
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
#include "randombytes_system.c"

static void print_fp2_hex(const char *label, const fp2_t *a)
{
    uint8_t buf[64];
    fp2_encode(buf, a);
    printf("%s = ", label);
    for (int i = 0; i < 64; i++)
        printf("%02x", buf[i]);
    printf("\n");
}

static void print_point_hex(const char *label, const ec_point_t *p)
{
    char lbl[128];
    snprintf(lbl, sizeof(lbl), "%s.x", label);
    print_fp2_hex(lbl, &p->x);
    snprintf(lbl, sizeof(lbl), "%s.z", label);
    print_fp2_hex(lbl, &p->z);
}

static void print_digit_array(const char *label, const digit_t *arr, int nwords)
{
    printf("%s = ", label);
    for (int w = 0; w < nwords; w++) {
        uint64_t val = arr[w];
        for (int b = 0; b < 8; b++) {
            printf("%02x", (unsigned)(val & 0xff));
            val >>= 8;
        }
    }
    printf("\n");
}

void
fp2_exp_2e(fp2_t *r, uint32_t e, const fp2_t *x)
{
    fp2_copy(r, x);
    for (uint32_t i = 0; i < e; i++) {
        fp2_sqr(r, r);
    }
}

int main(void)
{
    printf("=== Biextension Reference Test Trace (seed=0) ===\n\n");

    /* Init RNG with seed=0 */
    uint32_t seed[12] = {0};
    randombytes_init((unsigned char *)seed, NULL, 256);

    uint32_t e = TORSION_EVEN_POWER;
    ec_curve_t curve;
    ec_basis_t even_torsion;
    ec_point_t A24;

    /* E6: y^2 = x^3 + 6x^2 + x */
    ec_curve_init(&curve);
    fp2_set_small(&(curve.A), 6);
    fp2_set_one(&(curve.C));
    ec_curve_normalize_A24(&curve);
    copy_point(&A24, &curve.A24);

    /* Compute 2^e torsion on E6 */
    (void)ec_curve_to_basis_2f_to_hint(&even_torsion, &curve, e);

    printf("--- Basis points ---\n");
    print_point_hex("P", &even_torsion.P);
    print_point_hex("Q", &even_torsion.Q);
    print_point_hex("PmQ", &even_torsion.PmQ);

    /* Order checks */
    ec_point_t tmp, tmp2;
    ec_dbl_iter(&tmp, e, &even_torsion.P, &curve);
    assert(ec_is_zero(&tmp));
    ec_dbl_iter(&tmp, e, &even_torsion.Q, &curve);
    assert(ec_is_zero(&tmp));
    printf("order checks OK\n\n");

    /* Weil pairing */
    ec_point_t PQ;
    fp2_t weil_r, tate_r;
    xADD(&PQ, &even_torsion.P, &even_torsion.Q, &even_torsion.PmQ);
    weil(&weil_r, e, &even_torsion.P, &even_torsion.Q, &PQ, &curve);
    print_fp2_hex("weil", &weil_r);

    reduced_tate(&tate_r, e, &even_torsion.P, &even_torsion.Q, &PQ, &curve);
    print_fp2_hex("tate", &tate_r);

    /* Order of pairings */
    fp2_t one, r2;
    fp2_set_one(&one);
    fp2_exp_2e(&r2, e - 1, &weil_r);
    assert(!fp2_is_equal(&r2, &one));
    fp2_exp_2e(&r2, e, &weil_r);
    assert(fp2_is_equal(&r2, &one));
    printf("weil order OK\n");

    fp2_exp_2e(&r2, e - 1, &tate_r);
    assert(!fp2_is_equal(&r2, &one));
    fp2_exp_2e(&r2, e, &tate_r);
    assert(fp2_is_equal(&r2, &one));
    printf("tate order OK\n");

    /* Bilinearity */
    fp2_t weil_r2, weil_r3, rr1;
    weil(&weil_r2, e, &even_torsion.P, &even_torsion.Q, &even_torsion.PmQ, &curve);
    fp2_inv(&weil_r2);
    assert(fp2_is_equal(&weil_r, &weil_r2));
    printf("bilinearity OK\n\n");

    /* Double-bilinearity check */
    ec_point_t PP, QQ, PPQ, PQQ;
    xDBL_A24(&PP, &even_torsion.P, &A24, false);
    xDBL_A24(&QQ, &even_torsion.Q, &A24, false);
    xADD(&PPQ, &PQ, &even_torsion.P, &even_torsion.Q);
    xADD(&PQQ, &PQ, &even_torsion.Q, &even_torsion.P);
    weil(&weil_r2, e, &PP, &even_torsion.Q, &PPQ, &curve);
    weil(&weil_r3, e, &even_torsion.P, &QQ, &PQQ, &curve);
    assert(fp2_is_equal(&weil_r2, &weil_r3));
    fp2_sqr(&rr1, &weil_r);
    assert(fp2_is_equal(&rr1, &weil_r2));
    printf("double-bilinearity OK\n\n");

    /* ---- dlog tests ---- */
    printf("--- dlog tests ---\n");

    ec_basis_t BPQ, BRS;
    digit_t scal_r1[NWORDS_ORDER] = {0};
    digit_t scal_r2[NWORDS_ORDER] = {0};
    digit_t scal_s1[NWORDS_ORDER] = {0};
    digit_t scal_s2[NWORDS_ORDER] = {0};
    digit_t scal_d1[NWORDS_ORDER] = {0};
    digit_t scal_d2[NWORDS_ORDER] = {0};

    BPQ = even_torsion;
    BRS = even_torsion;

    randombytes((unsigned char *)scal_d1, (NWORDS_ORDER - 1) * sizeof(digit_t));
    randombytes((unsigned char *)scal_d2, (NWORDS_ORDER - 1) * sizeof(digit_t));
    randombytes((unsigned char *)scal_s1, (NWORDS_ORDER - 1) * sizeof(digit_t));
    randombytes((unsigned char *)scal_s2, (NWORDS_ORDER - 1) * sizeof(digit_t));

    print_digit_array("raw_d1", scal_d1, NWORDS_ORDER);
    print_digit_array("raw_d2", scal_d2, NWORDS_ORDER);
    print_digit_array("raw_s1", scal_s1, NWORDS_ORDER);
    print_digit_array("raw_s2", scal_s2, NWORDS_ORDER);

    scal_s1[0] = (scal_s1[0] & ((digit_t)(-1) - 1)) + 1;
    scal_d1[0] = (scal_d1[0] & ((digit_t)(-1) - 1));
    scal_s2[0] = (scal_s2[0] & ((digit_t)(-1) - 1)) + 1;
    scal_d2[0] = (scal_d2[0] & ((digit_t)(-1) - 1)) + 1;

    mp_add(scal_r1, scal_d1, scal_s1, NWORDS_ORDER);
    mp_add(scal_r2, scal_d2, scal_s2, NWORDS_ORDER);

    printf("\n--- After fixup ---\n");
    print_digit_array("scal_d1", scal_d1, NWORDS_ORDER);
    print_digit_array("scal_d2", scal_d2, NWORDS_ORDER);
    print_digit_array("scal_s1", scal_s1, NWORDS_ORDER);
    print_digit_array("scal_s2", scal_s2, NWORDS_ORDER);
    print_digit_array("scal_r1", scal_r1, NWORDS_ORDER);
    print_digit_array("scal_r2", scal_r2, NWORDS_ORDER);

    ec_biscalar_mul(&BRS.P, scal_r1, scal_r2, e, &BPQ, &curve);
    ec_biscalar_mul(&BRS.Q, scal_s1, scal_s2, e, &BPQ, &curve);
    ec_biscalar_mul(&BRS.PmQ, scal_d1, scal_d2, e, &BPQ, &curve);

    printf("\n--- BRS basis ---\n");
    print_point_hex("BRS_P", &BRS.P);
    print_point_hex("BRS_Q", &BRS.Q);
    print_point_hex("BRS_PmQ", &BRS.PmQ);

    /* Weil dlog */
    printf("\n--- Weil dlog ---\n");
    ec_dlog_2_weil(scal_r1, scal_r2, scal_s1, scal_s2, &BPQ, &BRS, &curve, e);
    print_digit_array("weil_r1", scal_r1, NWORDS_ORDER);
    print_digit_array("weil_r2", scal_r2, NWORDS_ORDER);
    print_digit_array("weil_s1", scal_s1, NWORDS_ORDER);
    print_digit_array("weil_s2", scal_s2, NWORDS_ORDER);

    /* Verify weil dlog */
    ec_biscalar_mul(&tmp, scal_r1, scal_r2, e, &BPQ, &curve);
    assert(ec_is_equal(&tmp, &BRS.P));
    ec_biscalar_mul(&tmp, scal_s1, scal_s2, e, &BPQ, &curve);
    assert(ec_is_equal(&tmp, &BRS.Q));
    printf("weil dlog verified OK\n");

    /* Tate dlog */
    printf("\n--- Tate dlog ---\n");
    ec_dlog_2_tate(scal_r1, scal_r2, scal_s1, scal_s2, &BPQ, &BRS, &curve, e);
    print_digit_array("tate_r1", scal_r1, NWORDS_ORDER);
    print_digit_array("tate_r2", scal_r2, NWORDS_ORDER);
    print_digit_array("tate_s1", scal_s1, NWORDS_ORDER);
    print_digit_array("tate_s2", scal_s2, NWORDS_ORDER);

    ec_biscalar_mul(&tmp, scal_r1, scal_r2, e, &BPQ, &curve);
    assert(ec_is_equal(&tmp, &BRS.P));
    ec_biscalar_mul(&tmp, scal_s1, scal_s2, e, &BPQ, &curve);
    assert(ec_is_equal(&tmp, &BRS.Q));
    printf("tate dlog verified OK\n");

    /* Partial torsion tate dlog (e=126) */
    printf("\n--- Tate partial dlog (e=126) ---\n");
    int e_full = TORSION_EVEN_POWER;
    int e_partial = 126;

    ec_basis_t BRS_partial;
    ec_dbl_iter(&BRS_partial.P, e_full - e_partial, &BRS.P, &curve);
    ec_dbl_iter(&BRS_partial.Q, e_full - e_partial, &BRS.Q, &curve);
    ec_dbl_iter(&BRS_partial.PmQ, e_full - e_partial, &BRS.PmQ, &curve);

    ec_dlog_2_tate(scal_r1, scal_r2, scal_s1, scal_s2, &BPQ, &BRS_partial, &curve, e_partial);
    print_digit_array("partial_r1", scal_r1, NWORDS_ORDER);
    print_digit_array("partial_r2", scal_r2, NWORDS_ORDER);
    print_digit_array("partial_s1", scal_s1, NWORDS_ORDER);
    print_digit_array("partial_s2", scal_s2, NWORDS_ORDER);

    ec_biscalar_mul(&tmp, scal_r1, scal_r2, e, &BPQ, &curve);
    ec_dbl_iter(&tmp, e_full - e_partial, &tmp, &curve);
    assert(ec_is_equal(&tmp, &BRS_partial.P));
    ec_biscalar_mul(&tmp, scal_s1, scal_s2, e, &BPQ, &curve);
    ec_dbl_iter(&tmp, e_full - e_partial, &tmp, &curve);
    assert(ec_is_equal(&tmp, &BRS_partial.Q));
    printf("partial tate dlog verified OK\n");

    printf("\nAll tests passed!\n");
    return 0;
}
