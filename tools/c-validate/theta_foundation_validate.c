/*
 * Cross-validation harness for sqisign-theta Groups 1, 2, 3.
 *
 * Validates: theta_precomputation, double_point, double_iter,
 *   is_product_theta_point, double_couple_point, hadamard,
 *   apply_isomorphism, base_change_matrix_multiplication.
 *
 * Build: tools/c-validate/build_theta_foundation.sh
 * Run:   tools/c-validate/theta_foundation_cval
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

static void print_point_hex(const char *label, const ec_point_t *p)
{
    char lbl[128];
    snprintf(lbl, sizeof(lbl), "%s.x", label);
    print_fp2_hex(lbl, &p->x);
    snprintf(lbl, sizeof(lbl), "%s.z", label);
    print_fp2_hex(lbl, &p->z);
}

int main(void)
{
    printf("=== Theta Foundation Cross-Validation ===\n\n");

    /* --- Section 1: Construct a theta structure from known values --- */
    printf("--- Section 1: Theta structure from known Fp2 values ---\n");

    /* Use small deterministic Fp2 values */
    theta_point_t null_pt;
    fp2_set_small(&null_pt.x, 3);
    fp2_set_small(&null_pt.y, 5);
    fp2_set_small(&null_pt.z, 7);
    fp2_set_small(&null_pt.t, 11);

    print_fp2_hex("null_x", &null_pt.x);
    print_fp2_hex("null_y", &null_pt.y);
    print_fp2_hex("null_z", &null_pt.z);
    print_fp2_hex("null_t", &null_pt.t);

    /* --- Section 2: Hadamard --- */
    printf("\n--- Section 2: Hadamard ---\n");
    theta_point_t had_out;
    hadamard(&had_out, &null_pt);
    print_fp2_hex("had_x", &had_out.x);
    print_fp2_hex("had_y", &had_out.y);
    print_fp2_hex("had_z", &had_out.z);
    print_fp2_hex("had_t", &had_out.t);

    /* Double hadamard = 4 * original (up to projective scaling) */
    theta_point_t had2;
    hadamard(&had2, &had_out);
    print_fp2_hex("had2_x", &had2.x);
    print_fp2_hex("had2_y", &had2.y);
    print_fp2_hex("had2_z", &had2.z);
    print_fp2_hex("had2_t", &had2.t);

    /* --- Section 3: to_squared_theta --- */
    printf("\n--- Section 3: to_squared_theta ---\n");
    theta_point_t sq_out;
    to_squared_theta(&sq_out, &null_pt);
    print_fp2_hex("sqth_x", &sq_out.x);
    print_fp2_hex("sqth_y", &sq_out.y);
    print_fp2_hex("sqth_z", &sq_out.z);
    print_fp2_hex("sqth_t", &sq_out.t);

    /* --- Section 4: theta_precomputation --- */
    printf("\n--- Section 4: theta_precomputation ---\n");
    theta_structure_t ts;
    memset(&ts, 0, sizeof(ts));
    ts.null_point = null_pt;
    ts.precomputation = false;
    theta_precomputation(&ts);

    print_fp2_hex("XYZ0", &ts.XYZ0);
    print_fp2_hex("YZT0", &ts.YZT0);
    print_fp2_hex("XZT0", &ts.XZT0);
    print_fp2_hex("XYT0", &ts.XYT0);
    print_fp2_hex("xyz0", &ts.xyz0);
    print_fp2_hex("yzt0", &ts.yzt0);
    print_fp2_hex("xzt0", &ts.xzt0);
    print_fp2_hex("xyt0", &ts.xyt0);

    /* --- Section 5: double_point --- */
    printf("\n--- Section 5: double_point ---\n");
    theta_point_t dbl_out;
    double_point(&dbl_out, &ts, &null_pt);
    print_fp2_hex("dbl_x", &dbl_out.x);
    print_fp2_hex("dbl_y", &dbl_out.y);
    print_fp2_hex("dbl_z", &dbl_out.z);
    print_fp2_hex("dbl_t", &dbl_out.t);

    /* --- Section 6: double_iter --- */
    printf("\n--- Section 6: double_iter (n=3) ---\n");
    theta_point_t iter_out;
    double_iter(&iter_out, &ts, &null_pt, 3);
    print_fp2_hex("iter3_x", &iter_out.x);
    print_fp2_hex("iter3_y", &iter_out.y);
    print_fp2_hex("iter3_z", &iter_out.z);
    print_fp2_hex("iter3_t", &iter_out.t);

    /* --- Section 7: is_product_theta_point --- */
    printf("\n--- Section 7: is_product_theta_point ---\n");
    /* (3,5,7,11): 3*11=33, 5*7=35, not equal => not product */
    uint32_t is_prod = is_product_theta_point(&null_pt);
    printf("is_product(3,5,7,11) = 0x%08x\n", is_prod);

    /* Construct a product point: (a, b, c, d) with a*d = b*c */
    /* e.g. (2, 3, 4, 6) since 2*6=12=3*4 */
    theta_point_t prod_pt;
    fp2_set_small(&prod_pt.x, 2);
    fp2_set_small(&prod_pt.y, 3);
    fp2_set_small(&prod_pt.z, 4);
    fp2_set_small(&prod_pt.t, 6);
    is_prod = is_product_theta_point(&prod_pt);
    printf("is_product(2,3,4,6) = 0x%08x\n", is_prod);

    /* --- Section 8: double_couple_point --- */
    printf("\n--- Section 8: double_couple_point ---\n");
    ec_curve_t E0;
    ec_curve_init(&E0);
    fp2_set_small(&E0.A, 6);
    fp2_set_one(&E0.C);
    ec_curve_normalize_A24(&E0);

    ec_basis_t basis;
    ec_curve_to_basis_2f_to_hint(&basis, &E0, TORSION_EVEN_POWER);

    theta_couple_curve_t E12;
    E12.E1 = E0;
    E12.E2 = E0;

    theta_couple_point_t cp;
    copy_point(&cp.P1, &basis.P);
    copy_point(&cp.P2, &basis.Q);

    theta_couple_point_t cp_dbl;
    double_couple_point(&cp_dbl, &cp, &E12);
    print_point_hex("cp_dbl_P1", &cp_dbl.P1);
    print_point_hex("cp_dbl_P2", &cp_dbl.P2);

    /* --- Section 9: double_couple_point_iter --- */
    printf("\n--- Section 9: double_couple_point_iter (n=3) ---\n");
    theta_couple_point_t cp_iter;
    double_couple_point_iter(&cp_iter, 3, &cp, &E12);
    print_point_hex("cp_iter_P1", &cp_iter.P1);
    print_point_hex("cp_iter_P2", &cp_iter.P2);

    /* --- Section 10: basis change matrix multiplication --- */
    printf("\n--- Section 10: basis change matrix multiplication ---\n");

    /* Construct M1 from SPLITTING_TRANSFORMS[0] and M2 from SPLITTING_TRANSFORMS[4] */
    basis_change_matrix_t M1, M2, M_prod;
    set_base_change_matrix_from_precomp(&M1, &SPLITTING_TRANSFORMS[0]);
    set_base_change_matrix_from_precomp(&M2, &SPLITTING_TRANSFORMS[4]);
    base_change_matrix_multiplication(&M_prod, &M1, &M2);

    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            char lbl[32];
            snprintf(lbl, sizeof(lbl), "M_prod[%d][%d]", i, j);
            print_fp2_hex(lbl, &M_prod.m[i][j]);
        }
    }

    /* --- Section 11: apply_isomorphism --- */
    printf("\n--- Section 11: apply_isomorphism ---\n");
    theta_point_t iso_out;
    apply_isomorphism(&iso_out, &M1, &null_pt);
    print_fp2_hex("iso_x", &iso_out.x);
    print_fp2_hex("iso_y", &iso_out.y);
    print_fp2_hex("iso_z", &iso_out.z);
    print_fp2_hex("iso_t", &iso_out.t);

    /* --- Section 12: apply_isomorphism_general with t=0 --- */
    printf("\n--- Section 12: apply_isomorphism_general (t=0) ---\n");
    theta_point_t pt_t0 = null_pt;
    fp2_set_zero(&pt_t0.t);
    theta_point_t iso_gen_out;
    apply_isomorphism_general(&iso_gen_out, &M1, &pt_t0, false);
    print_fp2_hex("isogen_x", &iso_gen_out.x);
    print_fp2_hex("isogen_y", &iso_gen_out.y);
    print_fp2_hex("isogen_z", &iso_gen_out.z);
    print_fp2_hex("isogen_t", &iso_gen_out.t);

    printf("\nAll sections done.\n");
    return 0;
}
