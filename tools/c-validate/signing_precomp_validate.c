/*
 * Cross-validation harness for signing-layer precomputed constants (Level 1).
 *
 * Prints every precomputed constant to stdout in a parseable format:
 *   - ibz_t values as decimal strings
 *   - fp limb arrays as comma-separated hex
 *   - ibz_mat_2x2_t as 4 decimal strings (row-major: [0][0],[0][1],[1][0],[1][1])
 *
 * Build: tools/c-validate/build_signing_precomp.sh
 * Run:   tools/c-validate/signing_precomp_cval > signing_precomp_expected.txt
 *
 * The Rust test crates/precomp/tests/c_crossvalidate_signing_precomp.rs
 * compares against this output.
 */

#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <stdbool.h>
#include <time.h>
#include <stdlib.h>

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

#include "quaternion.h"
#include "endomorphism_action.h"
#include "quaternion_data.h"
#include "torsion_constants.h"
#include "ec_params.h"
#include "e0_basis.h"
#include "quaternion_constants.h"

/* Print ibz_t as decimal string */
static void print_ibz(const char *label, const ibz_t *x) {
    /* Get the decimal string size */
    int sz = mpz_sizeinbase(*x, 10) + 2;
    char *buf = malloc(sz + 1);
    mpz_get_str(buf, 10, *x);
    printf("%s=%s\n", label, buf);
    free(buf);
}

/* Print fp limb array as comma-separated hex */
static void print_fp_limbs(const char *label, const digit_t *limbs, int nwords) {
    printf("%s=", label);
    for (int i = 0; i < nwords; i++) {
        if (i > 0) printf(",");
        printf("0x%llx", (unsigned long long)limbs[i]);
    }
    printf("\n");
}

/* Print fp2_t as two fp limb arrays */
static void print_fp2_limbs(const char *label, const fp2_t *a, int nwords) {
    char buf[256];
    snprintf(buf, sizeof(buf), "%s.re", label);
    print_fp_limbs(buf, a->re, nwords);
    snprintf(buf, sizeof(buf), "%s.im", label);
    print_fp_limbs(buf, a->im, nwords);
}

/* Print ibz_mat_2x2_t: 4 entries in row-major order */
static void print_mat2x2(const char *label, const ibz_mat_2x2_t *m) {
    char buf[256];
    snprintf(buf, sizeof(buf), "%s[0][0]", label);
    print_ibz(buf, &(*m)[0][0]);
    snprintf(buf, sizeof(buf), "%s[0][1]", label);
    print_ibz(buf, &(*m)[0][1]);
    snprintf(buf, sizeof(buf), "%s[1][0]", label);
    print_ibz(buf, &(*m)[1][0]);
    snprintf(buf, sizeof(buf), "%s[1][1]", label);
    print_ibz(buf, &(*m)[1][1]);
}

/* Print ibz_mat_4x4_t: 16 entries in row-major order */
static void print_mat4x4(const char *label, const ibz_mat_4x4_t *m) {
    char buf[256];
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            snprintf(buf, sizeof(buf), "%s[%d][%d]", label, i, j);
            print_ibz(buf, &(*m)[i][j]);
        }
    }
}

/* Print quat_alg_elem_t: denom + 4 coords */
static void print_elem(const char *label, const quat_alg_elem_t *e) {
    char buf[256];
    snprintf(buf, sizeof(buf), "%s.denom", label);
    print_ibz(buf, &e->denom);
    for (int i = 0; i < 4; i++) {
        snprintf(buf, sizeof(buf), "%s.coord[%d]", label, i);
        print_ibz(buf, &e->coord[i]);
    }
}

/* Print quat_lattice_t: denom + 4x4 basis */
static void print_lattice(const char *label, const quat_lattice_t *lat) {
    char buf[256];
    snprintf(buf, sizeof(buf), "%s.denom", label);
    print_ibz(buf, &lat->denom);
    snprintf(buf, sizeof(buf), "%s.basis", label);
    print_mat4x4(buf, &lat->basis);
}

int main(void)
{
    printf("=== Signing Precomp Cross-Validation (Level 1) ===\n\n");

    /* ---- Section 1: Torsion constants ---- */
    printf("# torsion_constants\n");
    printf("TORSION_2POWER_BYTES=32\n");
    print_ibz("TWO_TO_SECURITY_BITS", &TWO_TO_SECURITY_BITS);
    print_ibz("TORSION_PLUS_2POWER", &TORSION_PLUS_2POWER);
    print_ibz("SEC_DEGREE", &SEC_DEGREE);
    print_ibz("COM_DEGREE", &COM_DEGREE);
    printf("\n");

    /* ---- Section 2: EC params ---- */
    printf("# ec_params\n");
    printf("TORSION_EVEN_POWER=%d\n", TORSION_EVEN_POWER);
    printf("P_COFACTOR_FOR_2F_BITLENGTH=%d\n", P_COFACTOR_FOR_2F_BITLENGTH);
    printf("P_COFACTOR_FOR_2F=%llu\n", (unsigned long long)p_cofactor_for_2f[0]);
    printf("\n");

    /* ---- Section 3: Quaternion data ---- */
    printf("# quaternion_data\n");
    print_ibz("QUAT_prime_cofactor", &QUAT_prime_cofactor);
    print_ibz("QUATALG_P", &QUATALG_PINFTY.p);
    printf("\n");

    /* ---- Section 3a: Extremal orders ---- */
    for (int idx = 0; idx < 7; idx++) {
        char prefix[64];
        const quat_p_extremal_maximal_order_t *ord = &EXTREMAL_ORDERS[idx];

        snprintf(prefix, sizeof(prefix), "EXTREMAL_ORDER[%d]", idx);
        printf("# %s\n", prefix);

        printf("%s.q=%u\n", prefix, ord->q);

        char buf[128];
        snprintf(buf, sizeof(buf), "%s.order", prefix);
        print_lattice(buf, &ord->order);

        snprintf(buf, sizeof(buf), "%s.z", prefix);
        print_elem(buf, &ord->z);

        snprintf(buf, sizeof(buf), "%s.t", prefix);
        print_elem(buf, &ord->t);

        printf("\n");
    }

    /* ---- Section 3b: Connecting ideals ---- */
    for (int idx = 0; idx < 7; idx++) {
        char prefix[64];
        const quat_left_ideal_t *ideal = &CONNECTING_IDEALS[idx];

        snprintf(prefix, sizeof(prefix), "CONNECTING_IDEAL[%d]", idx);
        printf("# %s\n", prefix);

        char buf[128];
        snprintf(buf, sizeof(buf), "%s.norm", prefix);
        print_ibz(buf, &ideal->norm);
        snprintf(buf, sizeof(buf), "%s.lattice", prefix);
        print_lattice(buf, &ideal->lattice);

        printf("\n");
    }

    /* ---- Section 3c: Conjugating elements ---- */
    for (int idx = 0; idx < 7; idx++) {
        char prefix[64];
        snprintf(prefix, sizeof(prefix), "CONJUGATING_ELEM[%d]", idx);
        printf("# %s\n", prefix);
        print_elem(prefix, &CONJUGATING_ELEMENTS[idx]);
        printf("\n");
    }

    /* ---- Section 4: Endomorphism action ---- */
    for (int idx = 0; idx < 7; idx++) {
        char prefix[64];
        const curve_with_endomorphism_ring_t *entry = &CURVES_WITH_ENDOMORPHISMS[idx];

        snprintf(prefix, sizeof(prefix), "ENDOMORPHISM[%d]", idx);
        printf("# %s\n", prefix);

        /* Curve: A, C, A24 */
        char buf[128];
        snprintf(buf, sizeof(buf), "%s.curve.A", prefix);
        print_fp2_limbs(buf, &entry->curve.A, NWORDS_FIELD);
        snprintf(buf, sizeof(buf), "%s.curve.C", prefix);
        print_fp2_limbs(buf, &entry->curve.C, NWORDS_FIELD);
        snprintf(buf, sizeof(buf), "%s.curve.A24.x", prefix);
        print_fp2_limbs(buf, &entry->curve.A24.x, NWORDS_FIELD);
        snprintf(buf, sizeof(buf), "%s.curve.A24.z", prefix);
        print_fp2_limbs(buf, &entry->curve.A24.z, NWORDS_FIELD);

        /* Basis: P, Q, PmQ */
        snprintf(buf, sizeof(buf), "%s.basis.P.x", prefix);
        print_fp2_limbs(buf, &entry->basis_even.P.x, NWORDS_FIELD);
        snprintf(buf, sizeof(buf), "%s.basis.P.z", prefix);
        print_fp2_limbs(buf, &entry->basis_even.P.z, NWORDS_FIELD);
        snprintf(buf, sizeof(buf), "%s.basis.Q.x", prefix);
        print_fp2_limbs(buf, &entry->basis_even.Q.x, NWORDS_FIELD);
        snprintf(buf, sizeof(buf), "%s.basis.Q.z", prefix);
        print_fp2_limbs(buf, &entry->basis_even.Q.z, NWORDS_FIELD);
        snprintf(buf, sizeof(buf), "%s.basis.PmQ.x", prefix);
        print_fp2_limbs(buf, &entry->basis_even.PmQ.x, NWORDS_FIELD);
        snprintf(buf, sizeof(buf), "%s.basis.PmQ.z", prefix);
        print_fp2_limbs(buf, &entry->basis_even.PmQ.z, NWORDS_FIELD);

        /* Action matrices */
        snprintf(buf, sizeof(buf), "%s.action_i", prefix);
        print_mat2x2(buf, &entry->action_i);
        snprintf(buf, sizeof(buf), "%s.action_j", prefix);
        print_mat2x2(buf, &entry->action_j);
        snprintf(buf, sizeof(buf), "%s.action_k", prefix);
        print_mat2x2(buf, &entry->action_k);
        snprintf(buf, sizeof(buf), "%s.action_gen2", prefix);
        print_mat2x2(buf, &entry->action_gen2);
        snprintf(buf, sizeof(buf), "%s.action_gen3", prefix);
        print_mat2x2(buf, &entry->action_gen3);
        snprintf(buf, sizeof(buf), "%s.action_gen4", prefix);
        print_mat2x2(buf, &entry->action_gen4);

        printf("\n");
    }

    /* ---- Section 5: E0 basis ---- */
    printf("# e0_basis\n");
    print_fp2_limbs("BASIS_E0_PX", &BASIS_E0_PX, NWORDS_FIELD);
    print_fp2_limbs("BASIS_E0_QX", &BASIS_E0_QX, NWORDS_FIELD);
    printf("\n");

    /* ---- Section 6: Quaternion constants ---- */
    printf("# quaternion_constants\n");
    printf("QUAT_PRIMALITY_NUM_ITER=%d\n", QUAT_primality_num_iter);
    printf("QUAT_REPRES_BOUND_INPUT=%d\n", QUAT_repres_bound_input);
    printf("QUAT_EQUIV_BOUND_COEFF=%d\n", QUAT_equiv_bound_coeff);
    printf("FINDUV_BOX_SIZE=%d\n", FINDUV_box_size);
    printf("FINDUV_CUBE_SIZE=%d\n", FINDUV_cube_size);
    printf("\n");

    printf("# Done\n");
    return 0;
}
