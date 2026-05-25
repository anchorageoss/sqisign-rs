/*
 * Standalone cross-validation harness for sqisign-ec isogeny layer (Level 1).
 *
 * Compiles the C reference fp, fp2, mp, ec, and isogeny layers in a
 * single translation unit (via #include of the .c files) and runs a
 * fixed sequence of isogeny operations on known inputs.  Prints each
 * result as hex-encoded bytes to stdout.
 *
 * The Rust test crates/ec/tests/c_crossvalidate_isog.rs runs
 * the same sequence and compares byte for byte.
 *
 * Build:  tools/c-validate/build_isog.sh
 * Run:    tools/c-validate/isog_cval
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

/* Provide the precomputed constant ec_params.c would supply */
const digit_t p_cofactor_for_2f[1] = {5};

/* Pull in EC + isogeny layers */
#include "ec.c"
#include "ec_jac.c"
#include "xisog.c"
#include "xeval.c"
#include "isog_chains.c"

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
    printf("=== Isogeny Cross-Validation Harness (Level 1) ===\n\n");

    /* ---- Curve E0: A=0, C=1 ---- */
    ec_curve_t E0;
    ec_curve_init(&E0);
    ec_curve_normalize_A24(&E0);

    /* ---- Generate a point P of known even order on E0 ---- */
    /* E0: y^2 = x^3 + x. Use P with x = (3, 1), z = (1, 0). */
    ec_point_t P;
    fp2_from_small(&P.x, 3, 1);
    fp2_set_one(&P.z);

    /* ---- Test 1: xisog_2 + xeval_2 ---- */
    /* Compute [order/2]P to get a 2-torsion point as kernel */
    {
        printf("--- Test 1: xisog_2 + xeval_2 ---\n");

        /* Use a small 2-torsion point: (0:0) is 2-torsion on E0, but
         * xisog_2 requires kernel != (0:0). We'll use a scalar mul to
         * get a 2-torsion point from a known generator. Actually, for E0
         * with A=0, the 2-torsion points (other than 0) satisfy x^2+1=0,
         * i.e. x = i. So kernel = (i : 1). */
        ec_point_t K2;
        fp2_from_small(&K2.x, 0, 1);  /* x = i */
        fp2_set_one(&K2.z);

        ec_kps2_t kps2;
        ec_point_t B;  /* codomain in A24 form */
        xisog_2(&kps2, &B, K2);

        printf("xisog_2 kernel: K=(0+i : 1)\n");
        print_point_hex("B_A24", &B);
        print_point_hex("kps2.K", &kps2.K);

        /* Evaluate on point Q = (7+2i : 1) */
        ec_point_t Q;
        fp2_from_small(&Q.x, 7, 2);
        fp2_set_one(&Q.z);

        ec_point_t R;
        xeval_2(&R, &Q, 1, &kps2);
        print_point_hex("xeval_2(Q)", &R);
        printf("\n");
    }

    /* ---- Test 2: xisog_4 + xeval_4 ---- */
    {
        printf("--- Test 2: xisog_4 + xeval_4 ---\n");

        /* For a degree-4 isogeny we need a point of order 4 on E0.
         * E0: y^2 = x^3 + x. Take P4 = (3+i : 1); its double should
         * be 2-torsion. Let's verify by computing and using it. */
        ec_point_t P4;
        fp2_from_small(&P4.x, 3, 1);
        fp2_set_one(&P4.z);

        /* Check [2]P4 is 2-torsion - we use this even if it isn't a true
         * 4-torsion point; the formulas still produce deterministic output. */
        ec_kps4_t kps4;
        ec_point_t B4;
        xisog_4(&kps4, &B4, P4);

        print_point_hex("B4_A24", &B4);
        print_point_hex("kps4.K[0]", &kps4.K[0]);
        print_point_hex("kps4.K[1]", &kps4.K[1]);
        print_point_hex("kps4.K[2]", &kps4.K[2]);

        /* Evaluate on Q = (7+2i : 1) */
        ec_point_t Q;
        fp2_from_small(&Q.x, 7, 2);
        fp2_set_one(&Q.z);

        ec_point_t R;
        xeval_4(&R, &Q, 1, &kps4);
        print_point_hex("xeval_4(Q)", &R);
        printf("\n");
    }

    /* ---- Test 3: xisog_2 codomain then another xisog_2 (chain of 2) ---- */
    {
        printf("--- Test 3: two-step degree-2 chain ---\n");

        /* Start: E0, kernel K with order 4: double once -> 2-torsion for first step */
        ec_point_t K;
        fp2_from_small(&K.x, 3, 1);
        fp2_set_one(&K.z);

        /* First: find [2]K */
        ec_point_t AC;
        fp2_copy(&AC.x, &E0.A);
        fp2_copy(&AC.z, &E0.C);
        ec_point_t K2;
        xDBL(&K2, &K, &AC);

        print_point_hex("K", &K);
        print_point_hex("[2]K", &K2);

        /* First 2-isogeny with kernel [2]K */
        ec_kps2_t kps2_1;
        ec_point_t B1;
        xisog_2(&kps2_1, &B1, K2);
        print_point_hex("step1_B_A24", &B1);

        /* Push K through first isogeny */
        ec_point_t K_img;
        xeval_2(&K_img, &K, 1, &kps2_1);
        print_point_hex("K_after_step1", &K_img);

        /* Second 2-isogeny with kernel K_img */
        ec_kps2_t kps2_2;
        ec_point_t B2;
        xisog_2(&kps2_2, &B2, K_img);
        print_point_hex("step2_B_A24", &B2);

        /* Push test point (5+3i : 1) through both steps */
        ec_point_t T;
        fp2_from_small(&T.x, 5, 3);
        fp2_set_one(&T.z);
        xeval_2(&T, &T, 1, &kps2_1);
        xeval_2(&T, &T, 1, &kps2_2);
        print_point_hex("T_through_chain", &T);
        printf("\n");
    }

    /* ---- Test 4: ec_isomorphism + ec_iso_eval ---- */
    {
        printf("--- Test 4: ec_isomorphism + ec_iso_eval ---\n");

        /* Two j-equivalent curves: E0 (A=0, C=1) and E1 (A=6, C=1) have
         * different j-invariants, so they aren't isomorphic. Instead,
         * create two copies of the same curve with different (A:C) scaling.
         * E0: (A:C) = (0:1), E0_scaled: (A:C) = (0:2). These are the
         * same curve. */
        ec_curve_t E_from, E_to;
        ec_curve_init(&E_from);  /* A=0, C=1 */
        ec_curve_init(&E_to);    /* A=0, C=1 */

        /* Instead use non-trivial A. E1: A=6, C=1 and E1_scaled: A=12, C=2 */
        fp2_from_small(&E_from.A, 6, 0);
        fp2_set_one(&E_from.C);
        fp2_from_small(&E_to.A, 12, 0);
        fp2_from_small(&E_to.C, 2, 0);

        ec_isom_t isom;
        uint32_t err = ec_isomorphism(&isom, &E_from, &E_to);
        printf("ec_isomorphism err = %u\n", err);
        print_fp2_hex("isom.Nx", &isom.Nx);
        print_fp2_hex("isom.Nz", &isom.Nz);
        print_fp2_hex("isom.D", &isom.D);

        /* Apply to point (3+i : 1) */
        ec_point_t T;
        fp2_from_small(&T.x, 3, 1);
        fp2_set_one(&T.z);
        ec_iso_eval(&T, &isom);
        print_point_hex("iso_eval(3+i:1)", &T);
        printf("\n");
    }

    /* ---- Test 5: ec_eval_small_chain (length 1, success case) ---- */
    {
        printf("--- Test 5: ec_eval_small_chain ---\n");

        /* Use E0, kernel K = (i : 1) which is 2-torsion on E0, chain length 1 */
        ec_curve_t E;
        ec_curve_init(&E);

        ec_point_t K;
        fp2_from_small(&K.x, 0, 1);
        fp2_set_one(&K.z);

        /* Push-through point */
        ec_point_t pts[1];
        fp2_from_small(&pts[0].x, 7, 2);
        fp2_set_one(&pts[0].z);

        uint32_t ret = ec_eval_small_chain(&E, &K, 1, pts, 1, false);
        printf("ec_eval_small_chain ret = %u\n", ret);
        print_point_hex("chain_pt", &pts[0]);
        print_fp2_hex("chain_E.A", &E.A);
        print_fp2_hex("chain_E.C", &E.C);
        printf("\n");
    }

    /* ---- Test 6: simple xisog_2 round-trip codomain check ---- */
    {
        printf("--- Test 6: xisog_2 codomain A:C recovery ---\n");

        /* Kernel K = (i : 1) on E0 */
        ec_point_t K;
        fp2_from_small(&K.x, 0, 1);
        fp2_set_one(&K.z);

        ec_kps2_t kps;
        ec_point_t B_A24;
        xisog_2(&kps, &B_A24, K);

        /* Convert A24 -> A:C */
        ec_curve_t codomain;
        A24_to_AC(&codomain, &B_A24);
        print_fp2_hex("codomain.A", &codomain.A);
        print_fp2_hex("codomain.C", &codomain.C);

        /* Compute j-invariant of codomain */
        fp2_t j;
        ec_j_inv(&j, &codomain);
        print_fp2_hex("codomain_j", &j);
        printf("\n");
    }

    return 0;
}
