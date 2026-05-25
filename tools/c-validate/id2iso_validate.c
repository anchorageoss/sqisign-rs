/*
 * Cross-validation harness for id2iso mp_invert_matrix.
 *
 * This harness tests the digit-level matrix inversion used by
 * the verification path. It links only against mp.c (no quaternion deps).
 *
 * Output: hex-encoded results that the Rust test compares against.
 */

#include <stdio.h>
#include <string.h>
#include <inttypes.h>

#include <tutil.h>
#include <fp_constants.h>
#include <mp.h>

/* Print a digit array as hex (limb 0 first, each limb 16 hex chars). */
static void
print_digits(const digit_t *x, int nwords)
{
    for (int i = 0; i < nwords; i++) {
        printf("%016" PRIx64, x[i]);
    }
}

/* ===== Test 1: mp_invert_matrix, small values ===== */
static void
test_mp_invert_matrix_small(void)
{
    printf("=== mp_invert_matrix_small ===\n");

    /* Matrix [[3, 1], [2, 5]], det = 15-2 = 13 (odd, invertible mod 2^e) */
    digit_t r1[NWORDS_ORDER] = {3, 0, 0, 0};
    digit_t r2[NWORDS_ORDER] = {1, 0, 0, 0};
    digit_t s1[NWORDS_ORDER] = {2, 0, 0, 0};
    digit_t s2[NWORDS_ORDER] = {5, 0, 0, 0};
    int e = 256;

    printf("INPUT r1="); print_digits(r1, NWORDS_ORDER); printf("\n");
    printf("INPUT r2="); print_digits(r2, NWORDS_ORDER); printf("\n");
    printf("INPUT s1="); print_digits(s1, NWORDS_ORDER); printf("\n");
    printf("INPUT s2="); print_digits(s2, NWORDS_ORDER); printf("\n");

    mp_invert_matrix(r1, r2, s1, s2, e, NWORDS_ORDER);

    printf("OUTPUT r1="); print_digits(r1, NWORDS_ORDER); printf("\n");
    printf("OUTPUT r2="); print_digits(r2, NWORDS_ORDER); printf("\n");
    printf("OUTPUT s1="); print_digits(s1, NWORDS_ORDER); printf("\n");
    printf("OUTPUT s2="); print_digits(s2, NWORDS_ORDER); printf("\n");

    /* Verify: original * inverse should give identity mod 2^e */
    digit_t a[NWORDS_ORDER] = {3, 0, 0, 0};
    digit_t b[NWORDS_ORDER] = {1, 0, 0, 0};
    digit_t c[NWORDS_ORDER] = {2, 0, 0, 0};
    digit_t d[NWORDS_ORDER] = {5, 0, 0, 0};

    digit_t check00[NWORDS_ORDER], check01[NWORDS_ORDER];
    digit_t check10[NWORDS_ORDER], check11[NWORDS_ORDER];
    digit_t tmp1[NWORDS_ORDER], tmp2[NWORDS_ORDER];

    /* check00 = a*r1 + b*s1 mod 2^e */
    mp_mul(tmp1, a, r1, NWORDS_ORDER);
    mp_mul(tmp2, b, s1, NWORDS_ORDER);
    mp_add(check00, tmp1, tmp2, NWORDS_ORDER);
    mp_mod_2exp(check00, e, NWORDS_ORDER);

    /* check01 = a*r2 + b*s2 mod 2^e */
    mp_mul(tmp1, a, r2, NWORDS_ORDER);
    mp_mul(tmp2, b, s2, NWORDS_ORDER);
    mp_add(check01, tmp1, tmp2, NWORDS_ORDER);
    mp_mod_2exp(check01, e, NWORDS_ORDER);

    /* check10 = c*r1 + d*s1 mod 2^e */
    mp_mul(tmp1, c, r1, NWORDS_ORDER);
    mp_mul(tmp2, d, s1, NWORDS_ORDER);
    mp_add(check10, tmp1, tmp2, NWORDS_ORDER);
    mp_mod_2exp(check10, e, NWORDS_ORDER);

    /* check11 = c*r2 + d*s2 mod 2^e */
    mp_mul(tmp1, c, r2, NWORDS_ORDER);
    mp_mul(tmp2, d, s2, NWORDS_ORDER);
    mp_add(check11, tmp1, tmp2, NWORDS_ORDER);
    mp_mod_2exp(check11, e, NWORDS_ORDER);

    printf("CHECK [0][0]="); print_digits(check00, NWORDS_ORDER); printf("\n");
    printf("CHECK [0][1]="); print_digits(check01, NWORDS_ORDER); printf("\n");
    printf("CHECK [1][0]="); print_digits(check10, NWORDS_ORDER); printf("\n");
    printf("CHECK [1][1]="); print_digits(check11, NWORDS_ORDER); printf("\n");
}

/* ===== Test 2: mp_invert_matrix, larger values ===== */
static void
test_mp_invert_matrix_large(void)
{
    printf("\n=== mp_invert_matrix_large ===\n");

    digit_t r1[NWORDS_ORDER] = {0x123456789abcdef1ULL, 0x0fedcba987654321ULL, 0, 0};
    digit_t r2[NWORDS_ORDER] = {0x1111111111111111ULL, 0, 0, 0};
    digit_t s1[NWORDS_ORDER] = {0x2222222222222222ULL, 0, 0, 0};
    digit_t s2[NWORDS_ORDER] = {0x3333333333333333ULL, 0x4444444444444444ULL, 0, 0};
    int e = 256;

    printf("INPUT r1="); print_digits(r1, NWORDS_ORDER); printf("\n");
    printf("INPUT r2="); print_digits(r2, NWORDS_ORDER); printf("\n");
    printf("INPUT s1="); print_digits(s1, NWORDS_ORDER); printf("\n");
    printf("INPUT s2="); print_digits(s2, NWORDS_ORDER); printf("\n");

    mp_invert_matrix(r1, r2, s1, s2, e, NWORDS_ORDER);

    printf("OUTPUT r1="); print_digits(r1, NWORDS_ORDER); printf("\n");
    printf("OUTPUT r2="); print_digits(r2, NWORDS_ORDER); printf("\n");
    printf("OUTPUT s1="); print_digits(s1, NWORDS_ORDER); printf("\n");
    printf("OUTPUT s2="); print_digits(s2, NWORDS_ORDER); printf("\n");
}

/* ===== Test 3: mp_invert_matrix, e=128 (smaller modulus) ===== */
static void
test_mp_invert_matrix_128(void)
{
    printf("\n=== mp_invert_matrix_128 ===\n");

    digit_t r1[NWORDS_ORDER] = {7, 0, 0, 0};
    digit_t r2[NWORDS_ORDER] = {3, 0, 0, 0};
    digit_t s1[NWORDS_ORDER] = {4, 0, 0, 0};
    digit_t s2[NWORDS_ORDER] = {9, 0, 0, 0};
    int e = 128;

    printf("INPUT r1="); print_digits(r1, NWORDS_ORDER); printf("\n");
    printf("INPUT r2="); print_digits(r2, NWORDS_ORDER); printf("\n");
    printf("INPUT s1="); print_digits(s1, NWORDS_ORDER); printf("\n");
    printf("INPUT s2="); print_digits(s2, NWORDS_ORDER); printf("\n");

    mp_invert_matrix(r1, r2, s1, s2, e, NWORDS_ORDER);

    printf("OUTPUT r1="); print_digits(r1, NWORDS_ORDER); printf("\n");
    printf("OUTPUT r2="); print_digits(r2, NWORDS_ORDER); printf("\n");
    printf("OUTPUT s1="); print_digits(s1, NWORDS_ORDER); printf("\n");
    printf("OUTPUT s2="); print_digits(s2, NWORDS_ORDER); printf("\n");
}

int
main(void)
{
    test_mp_invert_matrix_small();
    test_mp_invert_matrix_large();
    test_mp_invert_matrix_128();
    return 0;
}
