/*
 * Standalone cross-validation harness for sqisign-precomp Level 1 constants.
 *
 * Prints each precomputed constant as hex-encoded bytes to stdout.
 * The Rust test crates/precomp/tests/c_crossvalidate_precomp.rs
 * compares byte for byte.
 *
 * Build:  tools/c-validate/build_precomp.sh
 * Run:    tools/c-validate/precomp_cval
 */

#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <stdbool.h>
#include <time.h>

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

/* Provide the precomputed constant ec_params.c would supply */
#include "ec_params.c"

/* Pull in EC layer for j-invariant computation */
#include "ec.c"
#include "ec_jac.c"

/* Pull in basis constants */
#include "e0_basis.c"

/* ---------- helpers ---------- */

static void print_hex(const char *label, const uint8_t *buf, int len)
{
    printf("%s = ", label);
    for (int i = 0; i < len; i++)
        printf("%02x", buf[i]);
    printf("\n");
}

static void print_fp2_hex(const char *label, const fp2_t *a)
{
    uint8_t buf[64];
    fp2_encode(buf, a);
    print_hex(label, buf, 64);
}

int main(void)
{
    printf("=== Precomp Cross-Validation Harness (Level 1) ===\n\n");

    /* ---- Scalar constants ---- */
    printf("TORSION_EVEN_POWER = %d\n", TORSION_EVEN_POWER);
    printf("P_COFACTOR_FOR_2F_BITLENGTH = %d\n", P_COFACTOR_FOR_2F_BITLENGTH);
    printf("p_cofactor_for_2f[0] = %lu\n", (unsigned long)p_cofactor_for_2f[0]);
    printf("\n");

    /* ---- E0 curve: A=0, C=1 ---- */
    {
        ec_curve_t E0;
        ec_curve_init(&E0);

        print_fp2_hex("E0.A", &E0.A);
        print_fp2_hex("E0.C", &E0.C);

        /* j-invariant */
        fp2_t j;
        ec_j_inv(&j, &E0);
        print_fp2_hex("j(E0)", &j);
        printf("\n");
    }

    /* ---- Basis points ---- */
    {
        print_fp2_hex("BASIS_E0_PX", &BASIS_E0_PX);
        print_fp2_hex("BASIS_E0_QX", &BASIS_E0_QX);
        printf("\n");

        /* Encode and re-decode round-trip check */
        uint8_t px_bytes[64], qx_bytes[64];
        fp2_encode(px_bytes, &BASIS_E0_PX);
        fp2_encode(qx_bytes, &BASIS_E0_QX);

        fp2_t px_rt, qx_rt;
        fp2_decode(&px_rt, px_bytes);
        fp2_decode(&qx_rt, qx_bytes);

        print_fp2_hex("BASIS_E0_PX_rt", &px_rt);
        print_fp2_hex("BASIS_E0_QX_rt", &qx_rt);
        printf("\n");

        /* Verify basis points are on E0: y^2 = x^3 + x */
        fp2_t y_sq, t;

        /* Check PX is on E0 */
        fp2_sqr(&y_sq, &BASIS_E0_PX);
        fp2_mul(&y_sq, &y_sq, &BASIS_E0_PX);  /* x^3 */
        fp2_add(&y_sq, &y_sq, &BASIS_E0_PX);  /* x^3 + x */
        int px_on_curve = fp2_is_square(&y_sq);
        printf("PX on E0: %d\n", px_on_curve);

        /* Check QX is on E0 */
        fp2_sqr(&y_sq, &BASIS_E0_QX);
        fp2_mul(&y_sq, &y_sq, &BASIS_E0_QX);
        fp2_add(&y_sq, &y_sq, &BASIS_E0_QX);
        int qx_on_curve = fp2_is_square(&y_sq);
        printf("QX on E0: %d\n", qx_on_curve);
        printf("\n");
    }

    /* ---- Cofactor multiplication test ----
     * Multiply a basis point by the cofactor to get a point in the 2^f
     * torsion subgroup. Then double f times to get the identity. */
    {
        ec_curve_t E0;
        ec_curve_init(&E0);
        ec_curve_normalize_A24(&E0);

        ec_point_t P;
        fp2_copy(&P.x, &BASIS_E0_PX);
        fp2_set_one(&P.z);

        /* Multiply by cofactor */
        ec_point_t cofP;
        xMUL(&cofP, &P, p_cofactor_for_2f, P_COFACTOR_FOR_2F_BITLENGTH, &E0);
        print_fp2_hex("cof*PX.x", &cofP.x);
        print_fp2_hex("cof*PX.z", &cofP.z);

        /* Double TORSION_EVEN_POWER times -> should be identity (z=0) */
        ec_point_t test = cofP;
        for (int i = 0; i < TORSION_EVEN_POWER; i++) {
            xDBL_A24(&test, &test, &E0.A24, true);
        }
        int is_identity = fp2_is_zero(&test.z);
        printf("cof*P doubled %d times is identity: %d\n", TORSION_EVEN_POWER, is_identity);
        printf("\n");
    }

    return 0;
}
