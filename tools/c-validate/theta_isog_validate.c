/*
 * Cross-validation harness for sqisign-theta Groups 4 and 5.
 *
 * Validates: gluing_change_of_basis, gluing_compute, gluing_eval_point,
 *   gluing_eval_point_special_case, theta_isogeny_compute,
 *   theta_isogeny_eval, verify_two_torsion.
 *
 * Build: tools/c-validate/build_theta_isog.sh
 * Run:   tools/c-validate/theta_isog_cval
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

/* Pull in biextension (needed for basis generation) */
#include "biextension.c"

/* Pull in HD layer */
#include "hd_splitting_transforms.c"
#include "theta_structure.c"
#include "theta_isogenies.c"
#include "hd.c"

static void print_fp2_hex(const char *label, const fp2_t *a)
{
    uint8_t buf[64];
    fp2_encode(buf, a);
    printf("%s = ", label);
    for (int i = 0; i < 64; i++)
        printf("%02x", buf[i]);
    printf("\n");
}

static void print_theta_point(const char *prefix, const theta_point_t *p)
{
    char lbl[128];
    snprintf(lbl, sizeof(lbl), "%s.x", prefix);
    print_fp2_hex(lbl, &p->x);
    snprintf(lbl, sizeof(lbl), "%s.y", prefix);
    print_fp2_hex(lbl, &p->y);
    snprintf(lbl, sizeof(lbl), "%s.z", prefix);
    print_fp2_hex(lbl, &p->z);
    snprintf(lbl, sizeof(lbl), "%s.t", prefix);
    print_fp2_hex(lbl, &p->t);
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

int main(void)
{
    printf("=== Theta Isogeny Cross-Validation ===\n\n");

    /*
     * Setup: E0 with A=6, C=1.
     * Generate full 2^TORSION_EVEN_POWER torsion basis.
     * Double down to 2^4 = 16 torsion (gives us 8-torsion for gluing
     * and the original 16-torsion maps to 2-torsion in theta model,
     * but we'll also generate separate test points).
     *
     * For the gluing, we use E1=E0, E2=E0 but with DIFFERENT kernel
     * generators on each curve (K1 on E1, K2 on E2) to avoid diagonal
     * degeneracy.
     */
    ec_curve_t E0;
    ec_curve_init(&E0);
    fp2_set_small(&E0.A, 6);
    fp2_set_one(&E0.C);
    ec_curve_normalize_A24(&E0);

    /* Generate full torsion basis */
    ec_basis_t full_basis;
    ec_curve_to_basis_2f_to_hint(&full_basis, &E0, TORSION_EVEN_POWER);

    /* Double down to 8-torsion */
    int dbl_count = TORSION_EVEN_POWER - 3;
    ec_point_t K1_8_xz = full_basis.P;
    ec_point_t K2_8_xz = full_basis.Q;
    ec_point_t K1m2_8_xz = full_basis.PmQ;
    for (int i = 0; i < dbl_count; i++) {
        xDBL_A24(&K1_8_xz, &K1_8_xz, &E0.A24, 0);
        xDBL_A24(&K2_8_xz, &K2_8_xz, &E0.A24, 0);
        xDBL_A24(&K1m2_8_xz, &K1m2_8_xz, &E0.A24, 0);
    }

    printf("--- Setup: 8-torsion points (XZ) ---\n");
    print_point_hex("K1_8", &K1_8_xz);
    print_point_hex("K2_8", &K2_8_xz);
    print_point_hex("K1m2_8", &K1m2_8_xz);

    /* Lift to Jacobian coordinates */
    ec_basis_t bas8 = { .P = K1_8_xz, .Q = K2_8_xz, .PmQ = K1m2_8_xz };
    jac_point_t xyK1, xyK2;
    uint32_t ok = lift_basis(&xyK1, &xyK2, &bas8, &E0);
    printf("lift_basis ok = %u\n", ok);

    printf("\n--- Setup: 8-torsion points (Jacobian) ---\n");
    print_jac_hex("xyK1", &xyK1);
    print_jac_hex("xyK2", &xyK2);

    /*
     * Use E1=E2=E0, but with CROSS kernel:
     *   xyK1_8 = (K1 on E1, K2 on E2)
     *   xyK2_8 = (K2 on E1, K1 on E2)
     * This avoids the diagonal degeneracy.
     */
    theta_couple_curve_t E12;
    E12.E1 = E0;
    E12.E2 = E0;

    theta_couple_jac_point_t xyK1_8, xyK2_8;
    xyK1_8.P1 = xyK1;
    xyK1_8.P2 = xyK2;
    xyK2_8.P1 = xyK2;
    xyK2_8.P2 = xyK1;

    /* --- Section 1: gluing_compute --- */
    printf("\n--- Section 1: gluing_compute ---\n");
    theta_gluing_t gluing;
    int ret = gluing_compute(&gluing, &E12, &xyK1_8, &xyK2_8, true);
    printf("gluing_compute ret = %d\n", ret);
    if (ret) {
        print_theta_point("gluing_codomain", &gluing.codomain);
        print_fp2_hex("gluing_imageK1_8.x", &gluing.imageK1_8.x);
        print_fp2_hex("gluing_imageK1_8.y", &gluing.imageK1_8.y);
        print_theta_point("gluing_precomp", &gluing.precomputation);

        /* Print the basis change matrix */
        printf("\n--- Section 1b: basis change matrix ---\n");
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                char lbl[32];
                snprintf(lbl, sizeof(lbl), "M[%d][%d]", i, j);
                print_fp2_hex(lbl, &gluing.M.m[i][j]);
            }
        }
    }

    /* --- Section 2: gluing_eval_point --- */
    printf("\n--- Section 2: gluing_eval_point ---\n");
    if (ret) {
        /*
         * Evaluate on the kernel generators themselves.
         * Phi(K1_8) should yield a point with structure (x:x:y:y).
         */
        theta_point_t eval_K1, eval_K2;
        gluing_eval_point(&eval_K1, &xyK1_8, &gluing);
        gluing_eval_point(&eval_K2, &xyK2_8, &gluing);
        print_theta_point("gluing_eval_K1", &eval_K1);
        print_theta_point("gluing_eval_K2", &eval_K2);
    }

    /* --- Section 3: gluing_eval_point_special_case --- */
    printf("\n--- Section 3: gluing_eval_point_special_case ---\n");
    if (ret) {
        theta_couple_point_t sc_pt;
        sc_pt.P1 = K1_8_xz;
        sc_pt.P2 = K2_8_xz;
        theta_point_t sc_img;
        int sc_ret = gluing_eval_point_special_case(&sc_img, &sc_pt, &gluing);
        printf("special_case ret = %d\n", sc_ret);
        if (sc_ret) {
            print_theta_point("special_case_img", &sc_img);
        }
    }

    /* --- Section 4: theta_isogeny_compute --- */
    printf("\n--- Section 4: theta_isogeny_compute ---\n");
    if (ret) {
        /*
         * For theta_isogeny_compute we need 8-torsion theta points
         * that are NOT in the kernel of the next isogeny.
         *
         * Use the chain approach: start with 2^6=64-torsion (halving
         * to get room), push through gluing, then double to 8-torsion.
         *
         * Actually simpler: use the full 2^TORSION_EVEN_POWER basis,
         * double down to 2^6 = 64, push through gluing (giving 2^3=8
         * torsion after the isogeny removes 3 bits). Wait no, the gluing
         * kernel is of degree 4, so post-gluing order = order/4 only for
         * points in the kernel. For points outside, order is preserved.
         *
         * Let's take a different approach: start with 2^5 = 32-torsion,
         * use the 8-torsion part for gluing, push the 32-torsion through
         * gluing to get 32-torsion theta points, then double twice to get
         * 8-torsion for theta_isogeny_compute.
         */

        /* Generate fresh 2^5 = 32-torsion from full basis */
        int dbl5 = TORSION_EVEN_POWER - 5;
        ec_point_t T1_32 = full_basis.P;
        ec_point_t T2_32 = full_basis.Q;
        ec_point_t T1m2_32 = full_basis.PmQ;
        for (int i = 0; i < dbl5; i++) {
            xDBL_A24(&T1_32, &T1_32, &E0.A24, 0);
            xDBL_A24(&T2_32, &T2_32, &E0.A24, 0);
            xDBL_A24(&T1m2_32, &T1m2_32, &E0.A24, 0);
        }

        /* Lift 32-torsion to Jacobian */
        ec_basis_t bas32 = { .P = T1_32, .Q = T2_32, .PmQ = T1m2_32 };
        jac_point_t xyT1_32, xyT2_32;
        lift_basis(&xyT1_32, &xyT2_32, &bas32, &E0);

        /* Same cross pattern for E1xE2 */
        theta_couple_jac_point_t xyT1_cp, xyT2_cp;
        xyT1_cp.P1 = xyT1_32;
        xyT1_cp.P2 = xyT2_32;
        xyT2_cp.P1 = xyT2_32;
        xyT2_cp.P2 = xyT1_32;

        /* Push through gluing */
        theta_point_t theta_T1, theta_T2;
        gluing_eval_point(&theta_T1, &xyT1_cp, &gluing);
        gluing_eval_point(&theta_T2, &xyT2_cp, &gluing);
        print_theta_point("theta_T1_32", &theta_T1);
        print_theta_point("theta_T2_32", &theta_T2);

        /* Double twice to get 8-torsion for theta_isogeny_compute */
        theta_structure_t ts;
        ts.null_point = gluing.codomain;
        ts.precomputation = false;
        theta_precomputation(&ts);

        theta_point_t T1_8_theta, T2_8_theta;
        double_iter(&T1_8_theta, &ts, &theta_T1, 2);
        double_iter(&T2_8_theta, &ts, &theta_T2, 2);
        print_theta_point("T1_8_theta", &T1_8_theta);
        print_theta_point("T2_8_theta", &T2_8_theta);

        theta_isogeny_t isog;
        int isog_ret = theta_isogeny_compute(&isog, &ts, &T1_8_theta,
                                              &T2_8_theta, false, false, false);
        printf("theta_isogeny_compute ret = %d\n", isog_ret);
        if (isog_ret) {
            print_theta_point("isog_codomain", &isog.codomain.null_point);
            print_theta_point("isog_precomp", &isog.precomputation);

            /* --- Section 5: theta_isogeny_eval --- */
            printf("\n--- Section 5: theta_isogeny_eval ---\n");
            theta_point_t eval_out;
            theta_isogeny_eval(&eval_out, &isog, &theta_T1);
            print_theta_point("isog_eval_out", &eval_out);
        }

        /* --- Section 6: theta_isogeny_compute_4 --- */
        printf("\n--- Section 6: theta_isogeny_compute_4 ---\n");
        {
            /* Get 4-torsion by doubling the 8-torsion once */
            theta_point_t T1_4, T2_4;
            double_point(&T1_4, &ts, &T1_8_theta);
            double_point(&T2_4, &ts, &T2_8_theta);
            print_theta_point("T1_4_theta", &T1_4);
            print_theta_point("T2_4_theta", &T2_4);

            theta_isogeny_t isog4;
            theta_isogeny_compute_4(&isog4, &ts, &T1_4, &T2_4, false, false);
            print_theta_point("isog4_codomain", &isog4.codomain.null_point);
            print_theta_point("isog4_precomp", &isog4.precomputation);

            /* Eval the _4 isogeny on the 32-torsion point */
            theta_point_t eval4_out;
            theta_isogeny_eval(&eval4_out, &isog4, &theta_T1);
            print_theta_point("isog4_eval_out", &eval4_out);
        }

        /* --- Section 7b: theta_isogeny_compute_2 --- */
        printf("\n--- Section 7b: theta_isogeny_compute_2 ---\n");
        {
            /* Get 2-torsion by doubling the 8-torsion twice */
            theta_point_t T1_2, T2_2;
            double_iter(&T1_2, &ts, &T1_8_theta, 2);
            double_iter(&T2_2, &ts, &T2_8_theta, 2);
            print_theta_point("T1_2_theta", &T1_2);
            print_theta_point("T2_2_theta", &T2_2);

            theta_isogeny_t isog2;
            theta_isogeny_compute_2(&isog2, &ts, &T1_2, &T2_2, false, false);
            print_theta_point("isog2_codomain", &isog2.codomain.null_point);
            print_theta_point("isog2_precomp", &isog2.precomputation);

            /* Eval the _2 isogeny */
            theta_point_t eval2_out;
            theta_isogeny_eval(&eval2_out, &isog2, &theta_T1);
            print_theta_point("isog2_eval_out", &eval2_out);
        }
    }

    /* --- Section 8: verify_two_torsion --- */
    printf("\n--- Section 8: verify_two_torsion ---\n");
    {
        /* Make 4-torsion and 2-torsion from 8-torsion */
        theta_couple_point_t K1_cp, K2_cp;
        K1_cp.P1 = K1_8_xz; K1_cp.P2 = K2_8_xz;
        K2_cp.P1 = K2_8_xz; K2_cp.P2 = K1_8_xz;

        theta_couple_point_t K1_4, K2_4, K1_2, K2_2;
        double_couple_point(&K1_4, &K1_cp, &E12);
        double_couple_point(&K2_4, &K2_cp, &E12);
        double_couple_point(&K1_2, &K1_4, &E12);
        double_couple_point(&K2_2, &K2_4, &E12);
        int v2t = verify_two_torsion(&K1_2, &K2_2, &E12);
        printf("verify_two_torsion(valid) = %d\n", v2t);

        /* Test with identical points (should fail) */
        int v2t_fail = verify_two_torsion(&K1_2, &K1_2, &E12);
        printf("verify_two_torsion(same) = %d\n", v2t_fail);
    }

    printf("\nAll sections done.\n");
    return 0;
}
