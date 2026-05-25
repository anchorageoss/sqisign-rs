/*
 * Standalone cross-validation harness for sqisign-ec biextension (Level 1).
 *
 * Exercises Weil/Tate pairings, pairing order checks, bilinearity,
 * discrete logs via Weil and Tate pairings, partial Tate dlog,
 * clear_cofac, and fp2_frob.
 *
 * Build:  tools/c-validate/build_biext.sh
 * Run:    tools/c-validate/biext_cval
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

/* ---------- helpers ---------- */

static void print_fp2_hex(const char *label, const fp2_t *a)
{
    uint8_t buf[64];
    fp2_encode(buf, a);
    printf("%s = ", label);
    for (int i = 0; i < 64; i++)
        printf("%02x", buf[i]);
    printf("\n");
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

static void fp2_from_small(fp2_t *out, int re, int im)
{
    fp2_set_zero(out);
    if (re >= 0) {
        fp_set_small(&out->re, (digit_t)re);
    } else {
        fp_t tmp;
        fp_set_small(&tmp, (digit_t)(-re));
        fp_neg(&out->re, &tmp);
    }
    if (im >= 0) {
        fp_set_small(&out->im, (digit_t)im);
    } else {
        fp_t tmp;
        fp_set_small(&tmp, (digit_t)(-im));
        fp_neg(&out->im, &tmp);
    }
}

int main(void)
{
    printf("=== Biextension Cross-Validation Harness (Level 1) ===\n\n");

    /* ---- Set up curve E1 (A=6) and its basis ---- */
    ec_curve_t E1;
    ec_curve_init(&E1);
    fp2_set_small(&E1.A, 6);
    fp2_set_one(&E1.C);
    ec_curve_normalize_A24(&E1);

    ec_basis_t basis;
    (void)ec_curve_to_basis_2f_to_hint(&basis, &E1, TORSION_EVEN_POWER);

    ec_point_t PQ;
    xADD(&PQ, &basis.P, &basis.Q, &basis.PmQ);

    /* ---- Test 1: Weil pairing on E1 (A=6) ---- */
    {
        printf("--- Test 1: Weil pairing on E1 (A=6) ---\n");
        fp2_t weil_result;
        ec_curve_t E1_copy;
        memcpy(&E1_copy, &E1, sizeof(ec_curve_t));
        weil(&weil_result, TORSION_EVEN_POWER, &basis.P, &basis.Q, &PQ, &E1_copy);
        print_fp2_hex("weil_E1", &weil_result);
        printf("\n");
    }

    /* Recompute weil_result for use in later tests (weil may modify curve) */
    fp2_t weil_result;
    {
        ec_curve_t E1_copy;
        memcpy(&E1_copy, &E1, sizeof(ec_curve_t));
        weil(&weil_result, TORSION_EVEN_POWER, &basis.P, &basis.Q, &PQ, &E1_copy);
    }

    /* ---- Test 2: Reduced Tate pairing on E1 (A=6) ---- */
    fp2_t tate_result;
    {
        printf("--- Test 2: Reduced Tate pairing on E1 (A=6) ---\n");
        ec_curve_t E1_copy;
        memcpy(&E1_copy, &E1, sizeof(ec_curve_t));
        reduced_tate(&tate_result, TORSION_EVEN_POWER, &basis.P, &basis.Q, &PQ, &E1_copy);
        print_fp2_hex("tate_E1", &tate_result);
        printf("\n");
    }

    /* ---- Test 3: Weil pairing order check ---- */
    {
        printf("--- Test 3: Weil pairing order check ---\n");
        fp2_t tmp;
        fp2_copy(&tmp, &weil_result);
        for (uint32_t i = 0; i < TORSION_EVEN_POWER; i++)
            fp2_sqr(&tmp, &tmp);
        print_fp2_hex("weil_pow_2e", &tmp);

        fp2_copy(&tmp, &weil_result);
        for (uint32_t i = 0; i < TORSION_EVEN_POWER - 1; i++)
            fp2_sqr(&tmp, &tmp);
        print_fp2_hex("weil_pow_2e_minus1", &tmp);
        printf("\n");
    }

    /* ---- Test 4: Tate pairing order check ---- */
    {
        printf("--- Test 4: Tate pairing order check ---\n");
        fp2_t tmp;
        fp2_copy(&tmp, &tate_result);
        for (uint32_t i = 0; i < TORSION_EVEN_POWER; i++)
            fp2_sqr(&tmp, &tmp);
        print_fp2_hex("tate_pow_2e", &tmp);

        fp2_copy(&tmp, &tate_result);
        for (uint32_t i = 0; i < TORSION_EVEN_POWER - 1; i++)
            fp2_sqr(&tmp, &tmp);
        print_fp2_hex("tate_pow_2e_minus1", &tmp);
        printf("\n");
    }

    /* ---- Test 5: Weil bilinearity check ---- */
    {
        printf("--- Test 5: Weil bilinearity ---\n");
        /* e(P, Q) computed with P-Q as difference (instead of P+Q)
         * gives the inverse pairing value. Verify:
         *   weil(P, Q, P-Q) = 1/weil(P, Q, P+Q) */
        fp2_t weil_neg;
        ec_curve_t E1_copy;
        memcpy(&E1_copy, &E1, sizeof(ec_curve_t));
        weil(&weil_neg, TORSION_EVEN_POWER, &basis.P, &basis.Q, &basis.PmQ, &E1_copy);
        fp2_inv(&weil_neg);
        print_fp2_hex("weil_neg_inv", &weil_neg);
        printf("\n");
    }

    /* ---- Set up BRS basis for dlog tests ---- */
    /* Use small known scalars: r1=7, r2=12, s1=3, s2=11
     * det = 7*11 - 12*3 = 77 - 36 = 41 (odd, invertible mod 2^248) */
    digit_t known_r1[NWORDS_ORDER] = {7, 0, 0, 0};
    digit_t known_r2[NWORDS_ORDER] = {12, 0, 0, 0};
    digit_t known_s1[NWORDS_ORDER] = {3, 0, 0, 0};
    digit_t known_s2[NWORDS_ORDER] = {11, 0, 0, 0};

    ec_basis_t BRS;
    {
        ec_curve_t E1_copy;
        memcpy(&E1_copy, &E1, sizeof(ec_curve_t));
        ec_biscalar_mul(&BRS.P, known_r1, known_r2, TORSION_EVEN_POWER, &basis, &E1_copy);
    }
    {
        ec_curve_t E1_copy;
        memcpy(&E1_copy, &E1, sizeof(ec_curve_t));
        ec_biscalar_mul(&BRS.Q, known_s1, known_s2, TORSION_EVEN_POWER, &basis, &E1_copy);
    }

    /* Compute R - S = [r1-s1]P + [r2-s2]Q = [4]P + [1]Q */
    digit_t diff_d1[NWORDS_ORDER] = {4, 0, 0, 0};
    digit_t diff_d2[NWORDS_ORDER] = {1, 0, 0, 0};
    {
        ec_curve_t E1_copy;
        memcpy(&E1_copy, &E1, sizeof(ec_curve_t));
        ec_biscalar_mul(&BRS.PmQ, diff_d1, diff_d2, TORSION_EVEN_POWER, &basis, &E1_copy);
    }

    /* ---- Test 6: Dlog with Weil pairing ---- */
    {
        printf("--- Test 6: Dlog with Weil pairing ---\n");
        digit_t rec_r1[NWORDS_ORDER], rec_r2[NWORDS_ORDER];
        digit_t rec_s1[NWORDS_ORDER], rec_s2[NWORDS_ORDER];

        /* ec_dlog_2_weil takes non-const PQ, so copy */
        ec_basis_t BPQ_copy;
        memcpy(&BPQ_copy, &basis, sizeof(ec_basis_t));
        ec_curve_t E1_copy;
        memcpy(&E1_copy, &E1, sizeof(ec_curve_t));

        ec_dlog_2_weil(rec_r1, rec_r2, rec_s1, rec_s2,
                       &BPQ_copy, &BRS, &E1_copy, TORSION_EVEN_POWER);

        print_digit_array("weil_dlog_r1", rec_r1, NWORDS_ORDER);
        print_digit_array("weil_dlog_r2", rec_r2, NWORDS_ORDER);
        print_digit_array("weil_dlog_s1", rec_s1, NWORDS_ORDER);
        print_digit_array("weil_dlog_s2", rec_s2, NWORDS_ORDER);
        printf("\n");
    }

    /* ---- Test 7: Dlog with Tate pairing ---- */
    {
        printf("--- Test 7: Dlog with Tate pairing ---\n");
        digit_t rec_r1[NWORDS_ORDER], rec_r2[NWORDS_ORDER];
        digit_t rec_s1[NWORDS_ORDER], rec_s2[NWORDS_ORDER];

        ec_curve_t E1_copy;
        memcpy(&E1_copy, &E1, sizeof(ec_curve_t));

        ec_dlog_2_tate(rec_r1, rec_r2, rec_s1, rec_s2,
                       &basis, &BRS, &E1_copy, TORSION_EVEN_POWER);

        print_digit_array("tate_dlog_r1", rec_r1, NWORDS_ORDER);
        print_digit_array("tate_dlog_r2", rec_r2, NWORDS_ORDER);
        print_digit_array("tate_dlog_s1", rec_s1, NWORDS_ORDER);
        print_digit_array("tate_dlog_s2", rec_s2, NWORDS_ORDER);
        printf("\n");
    }

    /* ---- Test 8: Tate partial dlog (e=126) ---- */
    {
        printf("--- Test 8: Tate partial dlog (e=126) ---\n");
        digit_t rec_r1[NWORDS_ORDER], rec_r2[NWORDS_ORDER];
        digit_t rec_s1[NWORDS_ORDER], rec_s2[NWORDS_ORDER];

        ec_basis_t BRS_partial;
        ec_curve_t E1_copy;
        memcpy(&E1_copy, &E1, sizeof(ec_curve_t));
        ec_dbl_iter(&BRS_partial.P, TORSION_EVEN_POWER - 126, &BRS.P, &E1_copy);

        memcpy(&E1_copy, &E1, sizeof(ec_curve_t));
        ec_dbl_iter(&BRS_partial.Q, TORSION_EVEN_POWER - 126, &BRS.Q, &E1_copy);

        memcpy(&E1_copy, &E1, sizeof(ec_curve_t));
        ec_dbl_iter(&BRS_partial.PmQ, TORSION_EVEN_POWER - 126, &BRS.PmQ, &E1_copy);

        memcpy(&E1_copy, &E1, sizeof(ec_curve_t));
        ec_dlog_2_tate(rec_r1, rec_r2, rec_s1, rec_s2,
                       &basis, &BRS_partial, &E1_copy, 126);

        print_digit_array("tate_partial_r1", rec_r1, NWORDS_ORDER);
        print_digit_array("tate_partial_r2", rec_r2, NWORDS_ORDER);
        print_digit_array("tate_partial_s1", rec_s1, NWORDS_ORDER);
        print_digit_array("tate_partial_s2", rec_s2, NWORDS_ORDER);
        printf("\n");
    }

    /* ---- Test 9: clear_cofac ---- */
    {
        printf("--- Test 9: clear_cofac ---\n");
        fp2_t cofac_input, cofac_output;
        fp2_from_small(&cofac_input, 3, 7);
        clear_cofac(&cofac_output, &cofac_input);
        print_fp2_hex("clear_cofac_3_7i", &cofac_output);
        printf("\n");
    }

    /* ---- Test 10: fp2_frob ---- */
    {
        printf("--- Test 10: fp2_frob ---\n");
        fp2_t frob_input, frob_output;
        fp2_from_small(&frob_input, 5, 11);
        fp2_frob(&frob_output, &frob_input);
        print_fp2_hex("fp2_frob_5_11i", &frob_output);
        printf("\n");
    }

    return 0;
}
