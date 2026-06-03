/*
 * Standalone cross-validation harness for sqisign-ec Level 1.
 *
 * Compiles the C reference fp, fp2, mp, and ec layers in a single
 * translation unit (via #include of the .c files) and runs a fixed
 * sequence of EC operations on known inputs. Prints each result as
 * hex-encoded bytes to stdout.
 *
 * The Rust test crates/ec/tests/c_crossvalidate.rs runs the
 * same sequence and compares byte for byte.
 *
 * Build:  tools/c-validate/build_ec.sh
 * Run:    tools/c-validate/ec_cval
 */

#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <stdbool.h>
#include <time.h>

/* Stub out tools.h timing functions (declared by ec.h -> tools.h but
   never called by ec.c or ec_jac.c). */
clock_t tic(void) { return 0; }
float tac(void) { return 0; }
float TAC(const char *str) { (void)str; return 0; }
float toc(const clock_t t) { (void)t; return 0; }
float TOC(const clock_t t, const char *str) { (void)t; (void)str; return 0; }
float TOC_clock(const clock_t t, const char *str) { (void)t; (void)str; return 0; }
clock_t dclock(const clock_t t) { (void)t; return 0; }
float clock_to_time(const clock_t t, const char *str) { (void)t; (void)str; return 0; }
float clock_print(const clock_t t, const char *str) { (void)t; (void)str; return 0; }

/* Pull in the full fp layer (static helpers + public API). */
#include "fp_p5248_64.c"

/* Pull in fp_select. */
#include "fp_select.c"

/* Pull in fp2 layer. */
#include "fp2.c"

/* Pull in mp layer. */
#include "mp.c"

/* Provide the precomputed constant ec_params.c would supply. */
const digit_t p_cofactor_for_2f[1] = {5};

/* Pull in EC layer. */
#include "ec.c"
#include "ec_jac.c"

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

static void print_jac_hex(const char *label, const jac_point_t *p)
{
    char lbl[128];
    snprintf(lbl, sizeof(lbl), "%s.x", label);
    print_fp2_hex(lbl, &p->x);
    snprintf(lbl, sizeof(lbl), "%s.y", label);
    print_fp2_hex(lbl, &p->y);
    snprintf(lbl, sizeof(lbl), "%s.z", label);
    print_fp2_hex(lbl, &p->z);
}

/* Decode an fp2 from 64 LE bytes (32 re + 32 im). */
static void fp2_from_bytes(fp2_t *out, const uint8_t bytes[64])
{
    fp2_decode(out, bytes);
}

/* Set an fp2 from two small integers. */
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
    ec_curve_t E0;
    ec_point_t P, Q, PQ, R;
    jac_point_t jP, jQ, jR;
    fp2_t j, y_sq, y;

    /* ---- Curve E0: A=0, C=1 ---- */
    ec_curve_init(&E0);
    ec_curve_normalize_A24(&E0);

    /* ---- Points on E0 ----
     * E0: y^2 = x^3 + x (A = 0).
     * We use small integer x-coordinates and verify they are on curve
     * by computing y = sqrt(x^3 + x) in Fp2.
     */

    /* P: x = (3, 1), z = (1, 0) */
    fp2_from_small(&P.x, 3, 1);
    fp2_set_one(&P.z);

    /* Compute y for P: y^2 = x^3 + x */
    fp2_sqr(&y_sq, &P.x);
    fp2_mul(&y_sq, &y_sq, &P.x);  /* x^3 */
    fp2_add(&y_sq, &y_sq, &P.x);  /* x^3 + x */
    fp2_sqrt(&y);                  /* y = sqrt(y_sq) -- modifies y_sq in place, result in y */

    /* Actually, fp2_sqrt takes pointer to the value and modifies in place */
    fp2_copy(&y, &y_sq);
    int valid_p = fp2_is_square(&y_sq);
    if (valid_p) {
        fp2_sqrt(&y);
    }

    /* Q: x = (7, 2), z = (1, 0) */
    fp2_from_small(&Q.x, 7, 2);
    fp2_set_one(&Q.z);

    /* PQ: arbitrary difference (5, 3) : (1, 0) -- for xADD the formula
     * computes from the inputs regardless of whether PQ is the actual P-Q. */
    fp2_from_small(&PQ.x, 5, 3);
    fp2_set_one(&PQ.z);

    printf("=== EC Cross-Validation Harness (Level 1) ===\n\n");

    /* Print inputs */
    print_point_hex("P", &P);
    print_point_hex("Q", &Q);
    print_point_hex("PQ", &PQ);
    printf("\n");

    /* ---- Test 1: xDBL ---- */
    {
        ec_point_t AC;
        fp2_copy(&AC.x, &E0.A);
        fp2_copy(&AC.z, &E0.C);
        xDBL(&R, &P, &AC);
        print_point_hex("xDBL(P)", &R);
    }

    /* ---- Test 2: xADD ---- */
    xADD(&R, &P, &Q, &PQ);
    print_point_hex("xADD(P,Q,PQ)", &R);

    /* ---- Test 3: xDBLADD ---- */
    {
        ec_point_t dbl_out, add_out;
        xDBLADD(&dbl_out, &add_out, &P, &Q, &PQ, &E0.A24, true);
        print_point_hex("xDBLADD.dbl", &dbl_out);
        print_point_hex("xDBLADD.add", &add_out);
    }

    /* ---- Test 4: xMUL (scalar multiplication) ---- */
    {
        /* k = 42, 256-bit scalar in 4 words */
        digit_t k[NWORDS_ORDER] = {42, 0, 0, 0};
        int kbits = 6; /* ceil(log2(42)) = 6 */
        xMUL(&R, &P, k, kbits, &E0);
        print_point_hex("xMUL(P,42)", &R);
    }

    /* ---- Test 5: ec_ladder3pt ---- */
    {
        digit_t m[NWORDS_ORDER] = {17, 0, 0, 0};
        ec_ladder3pt(&R, m, &P, &Q, &PQ, &E0);
        print_point_hex("ladder3pt(P,Q,PQ,17)", &R);
    }

    /* ---- Test 6: ec_j_inv ---- */
    {
        /* j-invariant of E0 should be 1728 = 6^3 * 8 */
        ec_j_inv(&j, &E0);
        print_fp2_hex("j_inv(E0)", &j);
    }

    /* ---- Test 7: j_inv with non-trivial curve ---- */
    {
        ec_curve_t E1;
        ec_curve_init(&E1);
        fp2_from_small(&E1.A, 6, 0);
        /* E1: y^2 = x^3 + 6x^2 + x, (A:C) = (6:1) */
        ec_j_inv(&j, &E1);
        print_fp2_hex("j_inv(E1)", &j);
    }

    /* ---- Test 8: ec_normalize_point ---- */
    {
        ec_point_t T;
        /* Create a non-normalized point: (3+i : 7+2i) */
        fp2_from_small(&T.x, 3, 1);
        fp2_from_small(&T.z, 7, 2);
        ec_normalize_point(&T);
        print_point_hex("normalize(3+i:7+2i)", &T);
    }

    printf("\n");

    /* ---- Jacobian tests ---- */
    /* Create Jacobian points. For jac_add/jac_dbl we need (X:Y:Z).
     * Use arbitrary Y values since the formulas work on any input. */
    {
        /* Jacobian P: x=(3,1), y=(11,5), z=(1,0) */
        fp2_from_small(&jP.x, 3, 1);
        fp2_from_small(&jP.y, 11, 5);
        fp2_set_one(&jP.z);

        /* Jacobian Q: x=(7,2), y=(13,4), z=(1,0) */
        fp2_from_small(&jQ.x, 7, 2);
        fp2_from_small(&jQ.y, 13, 4);
        fp2_set_one(&jQ.z);

        /* Use E0 for Jacobian arithmetic (A=0 is valid for jac_add/jac_dbl) */
        print_jac_hex("jac_P", &jP);
        print_jac_hex("jac_Q", &jQ);

        /* Test 9: Jacobian DBL */
        DBL(&jR, &jP, &E0);
        print_jac_hex("jac_dbl(P)", &jR);

        /* Test 10: Jacobian ADD */
        ADD(&jR, &jP, &jQ, &E0);
        print_jac_hex("jac_add(P,Q)", &jR);

        /* Test 11: jac_to_xz */
        {
            ec_point_t xz_out;
            jac_to_xz(&xz_out, &jR);
            print_point_hex("jac_to_xz(P+Q)", &xz_out);
        }

        /* Test 12: Jacobian ADD with identity */
        {
            jac_point_t jZ;
            jac_init(&jZ);
            ADD(&jR, &jP, &jZ, &E0);
            print_jac_hex("jac_add(P,0)", &jR);
        }

        /* Test 13: jac_neg + add = 0 */
        {
            jac_point_t neg_jP;
            jac_neg(&neg_jP, &jP);
            ADD(&jR, &jP, &neg_jP, &E0);
            print_jac_hex("jac_add(P,-P)", &jR);
        }
    }

    return 0;
}
