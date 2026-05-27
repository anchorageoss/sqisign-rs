/*
 * Cross-validation harness for sqisign-theta Groups 6 and 7.
 *
 * Validates: splitting_compute, theta_product_structure_to_elliptic_product,
 *   theta_point_to_montgomery_point.
 *
 * Build: tools/c-validate/build_theta_chain.sh
 * Run:   tools/c-validate/theta_chain_cval
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

#undef debug_print
#define debug_print(msg) ((void)0)

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

static void print_curve_hex(const char *label, const ec_curve_t *c)
{
    char lbl[128];
    snprintf(lbl, sizeof(lbl), "%s.A", label);
    print_fp2_hex(lbl, &c->A);
    snprintf(lbl, sizeof(lbl), "%s.C", label);
    print_fp2_hex(lbl, &c->C);
}

int main(void)
{
    printf("=== Theta Chain Cross-Validation ===\n\n");

    /*
     * Section 1: Construct a product theta structure and test
     * splitting_compute, theta_product_structure_to_elliptic_product,
     * and theta_point_to_montgomery_point.
     *
     * A product abelian surface E1 × E2 has a level-2 theta null point
     * of the form (a, b, c, 0) where one coordinate is zero, after
     * application of the correct splitting transform.
     *
     * We construct it by taking the codomain of a gluing isogeny from
     * E0 × E0 (with non-degenerate kernel), applying a series of theta
     * isogenies, and then testing splitting_compute on the result.
     *
     * Since E0 × E0 with a non-endomorphism kernel doesn't produce a
     * product, we instead construct the product theta structure directly
     * from the null point formula.
     *
     * For a product E1 × E2 with the identity splitting transform
     * (transform index 9 = identity), the null point is:
     *   null = (a, b, c, 0)
     * where a,b,c are derived from the fourth roots of theta functions.
     * We construct such a point and test splitting on it.
     */

    /* Simple product theta null point: (1, a, b, 0) in product form.
     * We use the formula from theta_product_structure_to_elliptic_product:
     *   E1: A1/C1 = -2(x^4+z^4)/(x^4-z^4)
     *   E2: A2/C2 = -2(x^4+y^4)/(x^4-y^4)
     *
     * Let's pick x = 1, y = 2, z = 3 (as Fp2 elements).
     * Then the null point is (1, 2, 3, 0), satisfying is_product.
     * E1: A1/C1 = -2(1+81)/(1-81) = -2*82/(-80) = 164/80 = 41/20
     * E2: A2/C2 = -2(1+16)/(1-16) = -2*17/(-15) = 34/15
     */
    printf("--- Section 1: splitting_compute on product theta structure ---\n");
    {
        theta_structure_t prod;
        fp2_set_small(&prod.null_point.x, 1);
        fp2_set_small(&prod.null_point.y, 2);
        fp2_set_small(&prod.null_point.z, 3);
        fp2_set_small(&prod.null_point.t, 6);
        prod.precomputation = 0;

        printf("is_product = %u\n", is_product_theta_point(&prod.null_point));

        theta_splitting_t split;
        bool ret = splitting_compute(&split, &prod, -1, false);
        printf("splitting_compute ret = %d\n", ret);

        if (ret) {
            printf("split null point:\n");
            print_theta_point("split_null", &split.B.null_point);
            printf("is_split_product = %u\n", is_product_theta_point(&split.B.null_point));
        }
    }

    /*
     * Section 2: theta_product_structure_to_elliptic_product.
     * Use the product null point (1, 2, 3, 0) directly.
     */
    printf("\n--- Section 2: theta_product_structure_to_elliptic_product ---\n");
    {
        theta_structure_t prod;
        fp2_set_small(&prod.null_point.x, 1);
        fp2_set_small(&prod.null_point.y, 2);
        fp2_set_small(&prod.null_point.z, 3);
        fp2_set_small(&prod.null_point.t, 6);
        prod.precomputation = 0;

        theta_couple_curve_t E12;
        int ret = theta_product_structure_to_elliptic_product(&E12, &prod);
        printf("product_to_elliptic ret = %d\n", ret);

        if (ret) {
            print_curve_hex("E1", &E12.E1);
            print_curve_hex("E2", &E12.E2);
        }
    }

    /*
     * Section 3: theta_point_to_montgomery_point.
     * Use the product null point (1, 2, 3, 0) as theta structure.
     * Construct a product theta point (a, b, c, 0) and convert.
     */
    printf("\n--- Section 3: theta_point_to_montgomery_point ---\n");
    {
        theta_structure_t prod;
        fp2_set_small(&prod.null_point.x, 1);
        fp2_set_small(&prod.null_point.y, 2);
        fp2_set_small(&prod.null_point.z, 3);
        fp2_set_small(&prod.null_point.t, 6);
        prod.precomputation = 0;

        /* A product point: (5, 7, 11, t) where t = 7*11/5.
         * 7*11 = 77, not divisible by 5. Use (5, 7, 10, 14) instead:
         * 5*14 = 70, 7*10 = 70 ✓ */
        theta_point_t P;
        fp2_set_small(&P.x, 5);
        fp2_set_small(&P.y, 7);
        fp2_set_small(&P.z, 10);
        fp2_set_small(&P.t, 14);

        theta_couple_point_t P12;
        int ret = theta_point_to_montgomery_point(&P12, &P, &prod);
        printf("point_to_montgomery ret = %d\n", ret);

        if (ret) {
            print_point_hex("P1", &P12.P1);
            print_point_hex("P2", &P12.P2);
        }

        /* Also test with a point where x=0, y=0 (fallback path).
         * Need xt = yz. With x=0, y=0: 0 = 0 ✓.
         * But theta_point_to_montgomery_point checks is_product first.
         * With (0, 0, z, t): xt=0, yz=0 ✓ as product point.
         * But x=0, y=0 triggers the fallback to use (z, t). */
        theta_point_t Q;
        fp2_set_zero(&Q.x);
        fp2_set_zero(&Q.y);
        fp2_set_small(&Q.z, 13);
        fp2_set_small(&Q.t, 17);

        theta_couple_point_t Q12;
        ret = theta_point_to_montgomery_point(&Q12, &Q, &prod);
        printf("point_to_montgomery (fallback) ret = %d\n", ret);

        if (ret) {
            print_point_hex("Q1_fb", &Q12.P1);
            print_point_hex("Q2_fb", &Q12.P2);
        }
    }

    /*
     * Section 4: splitting_compute with known zero_index.
     * The identity transform (index 9) should find the zero for our
     * product point (1, 2, 3, 0).
     */
    printf("\n--- Section 4: splitting_compute with zero_index ---\n");
    {
        theta_structure_t prod;
        fp2_set_small(&prod.null_point.x, 1);
        fp2_set_small(&prod.null_point.y, 2);
        fp2_set_small(&prod.null_point.z, 3);
        fp2_set_small(&prod.null_point.t, 6);
        prod.precomputation = 0;

        /* Test with zero_index = 9 (identity transform should match) */
        theta_splitting_t split;
        bool ret = splitting_compute(&split, &prod, 9, false);
        printf("splitting (zero_index=9) ret = %d\n", ret);

        /* Test with zero_index = 0 (should fail since index 0 doesn't match) */
        ret = splitting_compute(&split, &prod, 0, false);
        printf("splitting (zero_index=0) ret = %d\n", ret);
    }

    printf("\n=== Done ===\n");
    return 0;
}
