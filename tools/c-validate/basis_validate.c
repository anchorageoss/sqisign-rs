/*
 * Standalone cross-validation harness for sqisign-ec basis generation (Level 1).
 *
 * Build:  tools/c-validate/build_basis.sh
 * Run:    tools/c-validate/basis_cval
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

static void print_point_hex(const char *label, const ec_point_t *p)
{
    char lbl[128];
    snprintf(lbl, sizeof(lbl), "%s.x", label);
    print_fp2_hex(lbl, &p->x);
    snprintf(lbl, sizeof(lbl), "%s.z", label);
    print_fp2_hex(lbl, &p->z);
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
    printf("=== Basis Cross-Validation Harness (Level 1) ===\n\n");

    /* ---- Test 1: ec_basis_E0_2f for E0 (A=0) ---- */
    {
        printf("--- Test 1: ec_basis_E0_2f (A=0, f=248) ---\n");
        ec_curve_t E0;
        ec_curve_init(&E0);
        ec_curve_normalize_A24(&E0);

        ec_basis_t basis;
        uint8_t hint = ec_curve_to_basis_2f_to_hint(&basis, &E0, TORSION_EVEN_POWER);
        printf("hint = %u\n", hint);
        print_point_hex("P", &basis.P);
        print_point_hex("Q", &basis.Q);
        print_point_hex("PmQ", &basis.PmQ);
        printf("\n");
    }

    /* ---- Test 2: ec_basis_E0_2f for E0 with partial order (f=128) ---- */
    {
        printf("--- Test 2: ec_basis_E0_2f (A=0, f=128) ---\n");
        ec_curve_t E0;
        ec_curve_init(&E0);
        ec_curve_normalize_A24(&E0);

        ec_basis_t basis;
        uint8_t hint = ec_curve_to_basis_2f_to_hint(&basis, &E0, 128);
        printf("hint = %u\n", hint);
        print_point_hex("P128", &basis.P);
        print_point_hex("Q128", &basis.Q);
        print_point_hex("PmQ128", &basis.PmQ);
        printf("\n");
    }

    /* ---- Test 3: ec_recover_y ---- */
    {
        printf("--- Test 3: ec_recover_y ---\n");
        ec_curve_t E0;
        ec_curve_init(&E0);

        fp2_t px, y;
        fp2_from_small(&px, 3, 1);
        uint32_t ret = ec_recover_y(&y, &px, &E0);
        printf("ec_recover_y ret = %u\n", ret);
        print_fp2_hex("y(3+i)", &y);

        /* Also try on E1: A=6 */
        ec_curve_t E1;
        ec_curve_init(&E1);
        fp2_from_small(&E1.A, 6, 0);
        fp2_set_one(&E1.C);

        fp2_from_small(&px, 7, 2);
        ret = ec_recover_y(&y, &px, &E1);
        printf("ec_recover_y_E1 ret = %u\n", ret);
        print_fp2_hex("y_E1(7+2i)", &y);
        printf("\n");
    }

    /* ---- Test 4: ec_curve_to_basis_2f_to_hint on E1 (A=6) ---- */
    {
        printf("--- Test 4: basis_to_hint (A=6, f=248) ---\n");
        ec_curve_t E1;
        ec_curve_init(&E1);
        fp2_from_small(&E1.A, 6, 0);
        fp2_set_one(&E1.C);
        ec_curve_normalize_A24(&E1);

        ec_basis_t basis;
        uint8_t hint = ec_curve_to_basis_2f_to_hint(&basis, &E1, TORSION_EVEN_POWER);
        printf("hint = %u\n", hint);
        print_point_hex("P_E1", &basis.P);
        print_point_hex("Q_E1", &basis.Q);
        print_point_hex("PmQ_E1", &basis.PmQ);
        printf("\n");

        /* ---- Test 5: ec_curve_to_basis_2f_from_hint ---- */
        printf("--- Test 5: basis_from_hint (A=6, f=248, hint=%u) ---\n", hint);
        ec_curve_t E1b;
        ec_curve_init(&E1b);
        fp2_from_small(&E1b.A, 6, 0);
        fp2_set_one(&E1b.C);
        ec_curve_normalize_A24(&E1b);

        ec_basis_t basis2;
        int ok = ec_curve_to_basis_2f_from_hint(&basis2, &E1b, TORSION_EVEN_POWER, hint);
        printf("from_hint ok = %d\n", ok);
        print_point_hex("P_hint", &basis2.P);
        print_point_hex("Q_hint", &basis2.Q);
        print_point_hex("PmQ_hint", &basis2.PmQ);
        printf("\n");
    }

    /* ---- Test 6: is_on_curve ---- */
    {
        printf("--- Test 6: is_on_curve ---\n");
        ec_curve_t E0;
        ec_curve_init(&E0);

        /* (i : 1) should be on E0: y^2 = x^3 + x */
        fp2_t x;
        fp2_from_small(&x, 0, 1);

        /* E0 has C=1, so is_on_curve expects normalized curve */
        /* Actually is_on_curve is static in C... let's use ec_recover_y instead */
        fp2_t y;
        uint32_t ret = ec_recover_y(&y, &x, &E0);
        printf("(0+i) on E0: %u\n", ret);

        fp2_from_small(&x, 3, 1);
        ret = ec_recover_y(&y, &x, &E0);
        printf("(3+i) on E0: %u\n", ret);
        printf("\n");
    }

    /* ---- Test 7: lift_basis ---- */
    {
        printf("--- Test 7: lift_basis ---\n");
        ec_curve_t E1;
        ec_curve_init(&E1);
        fp2_from_small(&E1.A, 6, 0);
        fp2_set_one(&E1.C);
        ec_curve_normalize_A24(&E1);

        ec_basis_t basis;
        (void)ec_curve_to_basis_2f_to_hint(&basis, &E1, TORSION_EVEN_POWER);

        jac_point_t jP, jQ;
        uint32_t ret = lift_basis(&jP, &jQ, &basis, &E1);
        printf("lift_basis ret = %u\n", ret);

        /* Print Jacobian P */
        print_fp2_hex("jP.x", &jP.x);
        print_fp2_hex("jP.y", &jP.y);
        print_fp2_hex("jP.z", &jP.z);

        /* Print Jacobian Q */
        print_fp2_hex("jQ.x", &jQ.x);
        print_fp2_hex("jQ.y", &jQ.y);
        print_fp2_hex("jQ.z", &jQ.z);
        printf("\n");
    }

    return 0;
}
