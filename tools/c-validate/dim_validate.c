/*
 * Cross-validation harness for dim2.c, dim4.c, integers.c
 * → sqisign-quaternion::{dim2, dim4, integers}.
 *
 * Exercises vector/matrix/number-theory operations on hardcoded inputs
 * and prints hex-encoded results to stdout. The Rust test
 * (tests/c_crossvalidate_dim.rs) performs the identical operations
 * and compares byte-for-byte.
 *
 * Build: see build_dim.sh
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#include "intbig_internal.h"
#include "internal.h"
#include "quaternion.h"

static void print_ibz_hex(const char *label, const ibz_t *x) {
    int sz = ibz_size_in_base(x, 16) + 3;
    char buf[sz];
    memset(buf, 0, sz);
    ibz_convert_to_str(x, buf, 16);
    printf("%s=%s\n", label, buf);
}

static void print_vec2(const char *prefix, const ibz_vec_2_t *v) {
    char label[64];
    for (int i = 0; i < 2; i++) {
        snprintf(label, sizeof(label), "%s_%d", prefix, i);
        print_ibz_hex(label, &(*v)[i]);
    }
}

static void print_vec4(const char *prefix, const ibz_vec_4_t *v) {
    char label[64];
    for (int i = 0; i < 4; i++) {
        snprintf(label, sizeof(label), "%s_%d", prefix, i);
        print_ibz_hex(label, &(*v)[i]);
    }
}

static void print_mat2x2(const char *prefix, const ibz_mat_2x2_t *m) {
    char label[64];
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            snprintf(label, sizeof(label), "%s_%d%d", prefix, i, j);
            print_ibz_hex(label, &(*m)[i][j]);
        }
    }
}

static void print_mat4x4(const char *prefix, const ibz_mat_4x4_t *m) {
    char label[64];
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            snprintf(label, sizeof(label), "%s_%d%d", prefix, i, j);
            print_ibz_hex(label, &(*m)[i][j]);
        }
    }
}

int main(void) {
    /* ================================================================
     * Section 1: dim2, 2x2 matrix/vector arithmetic
     * ================================================================ */
    printf("# Section 1: dim2\n");

    /* 1a: ibz_vec_2_set */
    {
        ibz_vec_2_t v;
        ibz_vec_2_init(&v);
        ibz_vec_2_set(&v, 42, -17);
        print_vec2("vec2_set", &v);
        ibz_vec_2_finalize(&v);
    }

    /* 1b: ibz_mat_2x2_set and copy */
    {
        ibz_mat_2x2_t m, m2;
        ibz_mat_2x2_init(&m);
        ibz_mat_2x2_init(&m2);
        ibz_mat_2x2_set(&m, 3, -7, 11, 5);
        ibz_mat_2x2_copy(&m2, &m);
        print_mat2x2("mat2_set", &m2);
        ibz_mat_2x2_finalize(&m);
        ibz_mat_2x2_finalize(&m2);
    }

    /* 1c: ibz_mat_2x2_add */
    {
        ibz_mat_2x2_t a, b, sum;
        ibz_mat_2x2_init(&a);
        ibz_mat_2x2_init(&b);
        ibz_mat_2x2_init(&sum);
        ibz_mat_2x2_set(&a, 10, 20, 30, 40);
        ibz_mat_2x2_set(&b, 5, -3, 7, -9);
        ibz_mat_2x2_add(&sum, &a, &b);
        print_mat2x2("mat2_add", &sum);
        ibz_mat_2x2_finalize(&a);
        ibz_mat_2x2_finalize(&b);
        ibz_mat_2x2_finalize(&sum);
    }

    /* 1d: ibz_mat_2x2_det_from_ibz */
    {
        ibz_t a11, a12, a21, a22, det;
        ibz_init(&a11); ibz_init(&a12); ibz_init(&a21); ibz_init(&a22); ibz_init(&det);
        ibz_set(&a11, 3); ibz_set(&a12, 7);
        ibz_set(&a21, -2); ibz_set(&a22, 5);
        ibz_mat_2x2_det_from_ibz(&det, &a11, &a12, &a21, &a22);
        print_ibz_hex("det2x2", &det);
        ibz_finalize(&a11); ibz_finalize(&a12); ibz_finalize(&a21); ibz_finalize(&a22); ibz_finalize(&det);
    }

    /* 1e: ibz_mat_2x2_eval */
    {
        ibz_mat_2x2_t m;
        ibz_vec_2_t v, res;
        ibz_mat_2x2_init(&m);
        ibz_vec_2_init(&v);
        ibz_vec_2_init(&res);
        ibz_mat_2x2_set(&m, 3, -7, 11, 5);
        ibz_vec_2_set(&v, 4, -2);
        ibz_mat_2x2_eval(&res, &m, &v);
        print_vec2("mat2_eval", &res);
        ibz_mat_2x2_finalize(&m);
        ibz_vec_2_finalize(&v);
        ibz_vec_2_finalize(&res);
    }

    /* 1f: ibz_2x2_mul_mod */
    {
        ibz_mat_2x2_t a, b, prod;
        ibz_t mod;
        ibz_mat_2x2_init(&a);
        ibz_mat_2x2_init(&b);
        ibz_mat_2x2_init(&prod);
        ibz_init(&mod);
        ibz_mat_2x2_set(&a, 5, 3, -2, 7);
        ibz_mat_2x2_set(&b, 1, -4, 6, 2);
        ibz_set(&mod, 13);
        ibz_2x2_mul_mod(&prod, &a, &b, &mod);
        print_mat2x2("mat2_mulmod", &prod);
        ibz_mat_2x2_finalize(&a);
        ibz_mat_2x2_finalize(&b);
        ibz_mat_2x2_finalize(&prod);
        ibz_finalize(&mod);
    }

    /* 1g: ibz_mat_2x2_inv_mod (invertible) */
    {
        ibz_mat_2x2_t m, inv;
        ibz_t mod;
        ibz_mat_2x2_init(&m);
        ibz_mat_2x2_init(&inv);
        ibz_init(&mod);
        ibz_mat_2x2_set(&m, 5, 3, -2, 7);
        ibz_set(&mod, 13);
        int ok = ibz_mat_2x2_inv_mod(&inv, &m, &mod);
        printf("mat2_inv_ok=%d\n", ok);
        print_mat2x2("mat2_inv", &inv);
        ibz_mat_2x2_finalize(&m);
        ibz_mat_2x2_finalize(&inv);
        ibz_finalize(&mod);
    }

    /* 1h: ibz_mat_2x2_inv_mod (non-invertible: det=0 mod 7) */
    {
        ibz_mat_2x2_t m, inv;
        ibz_t mod;
        ibz_mat_2x2_init(&m);
        ibz_mat_2x2_init(&inv);
        ibz_init(&mod);
        ibz_mat_2x2_set(&m, 2, 3, 1, -2);
        ibz_set(&mod, 7);
        int ok = ibz_mat_2x2_inv_mod(&inv, &m, &mod);
        printf("mat2_inv_noninv_ok=%d\n", ok);
        ibz_mat_2x2_finalize(&m);
        ibz_mat_2x2_finalize(&inv);
        ibz_finalize(&mod);
    }

    /* ================================================================
     * Section 2: dim4, 4x4 matrix/vector arithmetic
     * ================================================================ */
    printf("# Section 2: dim4\n");

    /* 2a: ibz_vec_4_set */
    {
        ibz_vec_4_t v;
        ibz_vec_4_init(&v);
        ibz_vec_4_set(&v, 100, -200, 300, -400);
        print_vec4("vec4_set", &v);
        ibz_vec_4_finalize(&v);
    }

    /* 2b: ibz_vec_4_negate */
    {
        ibz_vec_4_t v, neg;
        ibz_vec_4_init(&v);
        ibz_vec_4_init(&neg);
        ibz_vec_4_set(&v, 1, -2, 3, -4);
        ibz_vec_4_negate(&neg, &v);
        print_vec4("vec4_neg", &neg);
        ibz_vec_4_finalize(&v);
        ibz_vec_4_finalize(&neg);
    }

    /* 2c: ibz_vec_4_add */
    {
        ibz_vec_4_t a, b, sum;
        ibz_vec_4_init(&a);
        ibz_vec_4_init(&b);
        ibz_vec_4_init(&sum);
        ibz_vec_4_set(&a, 10, 20, 30, 40);
        ibz_vec_4_set(&b, -5, 15, -25, 35);
        ibz_vec_4_add(&sum, &a, &b);
        print_vec4("vec4_add", &sum);
        ibz_vec_4_finalize(&a);
        ibz_vec_4_finalize(&b);
        ibz_vec_4_finalize(&sum);
    }

    /* 2d: ibz_vec_4_sub */
    {
        ibz_vec_4_t a, b, diff;
        ibz_vec_4_init(&a);
        ibz_vec_4_init(&b);
        ibz_vec_4_init(&diff);
        ibz_vec_4_set(&a, 10, 20, 30, 40);
        ibz_vec_4_set(&b, -5, 15, -25, 35);
        ibz_vec_4_sub(&diff, &a, &b);
        print_vec4("vec4_sub", &diff);
        ibz_vec_4_finalize(&a);
        ibz_vec_4_finalize(&b);
        ibz_vec_4_finalize(&diff);
    }

    /* 2e: ibz_vec_4_scalar_mul */
    {
        ibz_vec_4_t v, prod;
        ibz_t scalar;
        ibz_vec_4_init(&v);
        ibz_vec_4_init(&prod);
        ibz_init(&scalar);
        ibz_vec_4_set(&v, 3, -7, 11, -13);
        ibz_set(&scalar, 5);
        ibz_vec_4_scalar_mul(&prod, &scalar, &v);
        print_vec4("vec4_smul", &prod);
        ibz_vec_4_finalize(&v);
        ibz_vec_4_finalize(&prod);
        ibz_finalize(&scalar);
    }

    /* 2f: ibz_vec_4_scalar_div */
    {
        ibz_vec_4_t v, quot;
        ibz_t scalar;
        ibz_vec_4_init(&v);
        ibz_vec_4_init(&quot);
        ibz_init(&scalar);
        ibz_vec_4_set(&v, 15, -35, 55, -65);
        ibz_set(&scalar, 5);
        int ok = ibz_vec_4_scalar_div(&quot, &scalar, &v);
        printf("vec4_sdiv_ok=%d\n", ok);
        print_vec4("vec4_sdiv", &quot);
        ibz_vec_4_finalize(&v);
        ibz_vec_4_finalize(&quot);
        ibz_finalize(&scalar);
    }

    /* 2g: ibz_vec_4_content */
    {
        ibz_vec_4_t v;
        ibz_t content;
        ibz_vec_4_init(&v);
        ibz_init(&content);
        ibz_vec_4_set(&v, 12, -18, 24, -30);
        ibz_vec_4_content(&content, &v);
        print_ibz_hex("vec4_content", &content);
        ibz_vec_4_finalize(&v);
        ibz_finalize(&content);
    }

    /* 2h: ibz_vec_4_linear_combination */
    {
        ibz_vec_4_t a, b, lc;
        ibz_t ca, cb;
        ibz_vec_4_init(&a);
        ibz_vec_4_init(&b);
        ibz_vec_4_init(&lc);
        ibz_init(&ca);
        ibz_init(&cb);
        ibz_vec_4_set(&a, 1, 2, 3, 4);
        ibz_vec_4_set(&b, 5, 6, 7, 8);
        ibz_set(&ca, 3);
        ibz_set(&cb, -2);
        ibz_vec_4_linear_combination(&lc, &ca, &a, &cb, &b);
        print_vec4("vec4_lincomb", &lc);
        ibz_vec_4_finalize(&a);
        ibz_vec_4_finalize(&b);
        ibz_vec_4_finalize(&lc);
        ibz_finalize(&ca);
        ibz_finalize(&cb);
    }

    /* 2i: ibz_vec_4_is_zero */
    {
        ibz_vec_4_t v;
        ibz_vec_4_init(&v);
        ibz_vec_4_set(&v, 0, 0, 0, 0);
        printf("vec4_iszero_yes=%d\n", ibz_vec_4_is_zero(&v));
        ibz_vec_4_set(&v, 0, 0, 1, 0);
        printf("vec4_iszero_no=%d\n", ibz_vec_4_is_zero(&v));
        ibz_vec_4_finalize(&v);
    }

    /* 2j: ibz_mat_4x4_identity, is_identity, zero */
    {
        ibz_mat_4x4_t id, z;
        ibz_mat_4x4_init(&id);
        ibz_mat_4x4_init(&z);
        ibz_mat_4x4_identity(&id);
        printf("mat4_isid_yes=%d\n", ibz_mat_4x4_is_identity(&id));
        ibz_mat_4x4_zero(&z);
        printf("mat4_isid_no=%d\n", ibz_mat_4x4_is_identity(&z));
        ibz_mat_4x4_finalize(&id);
        ibz_mat_4x4_finalize(&z);
    }

    /* 2k: ibz_mat_4x4_mul */
    {
        ibz_mat_4x4_t a, b, prod;
        ibz_mat_4x4_init(&a);
        ibz_mat_4x4_init(&b);
        ibz_mat_4x4_init(&prod);

        /* a = [[1,2,3,4],[5,6,7,8],[9,10,11,12],[13,14,15,16]] */
        for (int i = 0; i < 4; i++)
            for (int j = 0; j < 4; j++)
                ibz_set(&a[i][j], i * 4 + j + 1);

        /* b = [[2,0,1,0],[0,3,0,1],[1,0,2,0],[0,1,0,3]] */
        ibz_mat_4x4_zero(&b);
        ibz_set(&b[0][0], 2); ibz_set(&b[0][2], 1);
        ibz_set(&b[1][1], 3); ibz_set(&b[1][3], 1);
        ibz_set(&b[2][0], 1); ibz_set(&b[2][2], 2);
        ibz_set(&b[3][1], 1); ibz_set(&b[3][3], 3);

        ibz_mat_4x4_mul(&prod, &a, &b);
        print_mat4x4("mat4_mul", &prod);

        ibz_mat_4x4_finalize(&a);
        ibz_mat_4x4_finalize(&b);
        ibz_mat_4x4_finalize(&prod);
    }

    /* 2l: ibz_mat_4x4_transpose */
    {
        ibz_mat_4x4_t m, mt;
        ibz_mat_4x4_init(&m);
        ibz_mat_4x4_init(&mt);
        for (int i = 0; i < 4; i++)
            for (int j = 0; j < 4; j++)
                ibz_set(&m[i][j], i * 4 + j + 1);
        ibz_mat_4x4_transpose(&mt, &m);
        print_mat4x4("mat4_trans", &mt);
        ibz_mat_4x4_finalize(&m);
        ibz_mat_4x4_finalize(&mt);
    }

    /* 2m: ibz_mat_4x4_scalar_mul */
    {
        ibz_mat_4x4_t m, prod;
        ibz_t scalar;
        ibz_mat_4x4_init(&m);
        ibz_mat_4x4_init(&prod);
        ibz_init(&scalar);
        for (int i = 0; i < 4; i++)
            for (int j = 0; j < 4; j++)
                ibz_set(&m[i][j], i * 4 + j + 1);
        ibz_set(&scalar, -3);
        ibz_mat_4x4_scalar_mul(&prod, &scalar, &m);
        print_mat4x4("mat4_smul", &prod);
        ibz_mat_4x4_finalize(&m);
        ibz_mat_4x4_finalize(&prod);
        ibz_finalize(&scalar);
    }

    /* 2n: ibz_mat_4x4_scalar_div */
    {
        ibz_mat_4x4_t m, quot;
        ibz_t scalar;
        ibz_mat_4x4_init(&m);
        ibz_mat_4x4_init(&quot);
        ibz_init(&scalar);
        for (int i = 0; i < 4; i++)
            for (int j = 0; j < 4; j++)
                ibz_set(&m[i][j], (i * 4 + j + 1) * 6);
        ibz_set(&scalar, 3);
        int ok = ibz_mat_4x4_scalar_div(&quot, &scalar, &m);
        printf("mat4_sdiv_ok=%d\n", ok);
        print_mat4x4("mat4_sdiv", &quot);
        ibz_mat_4x4_finalize(&m);
        ibz_mat_4x4_finalize(&quot);
        ibz_finalize(&scalar);
    }

    /* 2o: ibz_mat_4x4_negate */
    {
        ibz_mat_4x4_t m, neg;
        ibz_mat_4x4_init(&m);
        ibz_mat_4x4_init(&neg);
        for (int i = 0; i < 4; i++)
            for (int j = 0; j < 4; j++)
                ibz_set(&m[i][j], i * 4 + j + 1);
        ibz_mat_4x4_negate(&neg, &m);
        print_mat4x4("mat4_neg", &neg);
        ibz_mat_4x4_finalize(&m);
        ibz_mat_4x4_finalize(&neg);
    }

    /* 2p: ibz_mat_4x4_gcd */
    {
        ibz_mat_4x4_t m;
        ibz_t g;
        ibz_mat_4x4_init(&m);
        ibz_init(&g);
        for (int i = 0; i < 4; i++)
            for (int j = 0; j < 4; j++)
                ibz_set(&m[i][j], (i * 4 + j + 1) * 6);
        ibz_mat_4x4_gcd(&g, &m);
        print_ibz_hex("mat4_gcd", &g);
        ibz_mat_4x4_finalize(&m);
        ibz_finalize(&g);
    }

    /* 2q: ibz_mat_4x4_eval */
    {
        ibz_mat_4x4_t m;
        ibz_vec_4_t v, res;
        ibz_mat_4x4_init(&m);
        ibz_vec_4_init(&v);
        ibz_vec_4_init(&res);
        for (int i = 0; i < 4; i++)
            for (int j = 0; j < 4; j++)
                ibz_set(&m[i][j], i * 4 + j + 1);
        ibz_vec_4_set(&v, 1, -1, 2, -2);
        ibz_mat_4x4_eval(&res, &m, &v);
        print_vec4("mat4_eval", &res);
        ibz_mat_4x4_finalize(&m);
        ibz_vec_4_finalize(&v);
        ibz_vec_4_finalize(&res);
    }

    /* 2r: ibz_mat_4x4_eval_t */
    {
        ibz_mat_4x4_t m;
        ibz_vec_4_t v, res;
        ibz_mat_4x4_init(&m);
        ibz_vec_4_init(&v);
        ibz_vec_4_init(&res);
        for (int i = 0; i < 4; i++)
            for (int j = 0; j < 4; j++)
                ibz_set(&m[i][j], i * 4 + j + 1);
        ibz_vec_4_set(&v, 1, -1, 2, -2);
        ibz_mat_4x4_eval_t(&res, &v, &m);
        print_vec4("mat4_evalt", &res);
        ibz_mat_4x4_finalize(&m);
        ibz_vec_4_finalize(&v);
        ibz_vec_4_finalize(&res);
    }

    /* 2s: ibz_mat_4x4_equal */
    {
        ibz_mat_4x4_t a, b;
        ibz_mat_4x4_init(&a);
        ibz_mat_4x4_init(&b);
        for (int i = 0; i < 4; i++)
            for (int j = 0; j < 4; j++) {
                ibz_set(&a[i][j], i * 4 + j + 1);
                ibz_set(&b[i][j], i * 4 + j + 1);
            }
        printf("mat4_equal_yes=%d\n", ibz_mat_4x4_equal(&a, &b));
        ibz_set(&b[2][3], 999);
        printf("mat4_equal_no=%d\n", ibz_mat_4x4_equal(&a, &b));
        ibz_mat_4x4_finalize(&a);
        ibz_mat_4x4_finalize(&b);
    }

    /* 2t: ibz_mat_4x4_inv_with_det_as_denom */
    {
        ibz_mat_4x4_t m, inv;
        ibz_t det;
        ibz_mat_4x4_init(&m);
        ibz_mat_4x4_init(&inv);
        ibz_init(&det);

        /* Upper triangular: [[2,1,3,0],[0,4,0,1],[0,0,3,2],[0,0,0,2]], det=48 */
        ibz_mat_4x4_zero(&m);
        ibz_set(&m[0][0], 2); ibz_set(&m[0][1], 1); ibz_set(&m[0][2], 3);
        ibz_set(&m[1][1], 4); ibz_set(&m[1][3], 1);
        ibz_set(&m[2][2], 3); ibz_set(&m[2][3], 2);
        ibz_set(&m[3][3], 2);

        int ok = ibz_mat_4x4_inv_with_det_as_denom(&inv, &det, &m);
        printf("mat4_inv_ok=%d\n", ok);
        print_ibz_hex("mat4_inv_det", &det);
        print_mat4x4("mat4_inv", &inv);

        ibz_mat_4x4_finalize(&m);
        ibz_mat_4x4_finalize(&inv);
        ibz_finalize(&det);
    }

    /* 2u: quat_qf_eval */
    {
        ibz_mat_4x4_t qf;
        ibz_vec_4_t coord;
        ibz_t result;
        ibz_mat_4x4_init(&qf);
        ibz_vec_4_init(&coord);
        ibz_init(&result);

        /* qf = diag(1,1,3,3) representing x0^2+x1^2+3*x2^2+3*x3^2 */
        ibz_mat_4x4_zero(&qf);
        ibz_set(&qf[0][0], 1);
        ibz_set(&qf[1][1], 1);
        ibz_set(&qf[2][2], 3);
        ibz_set(&qf[3][3], 3);

        ibz_vec_4_set(&coord, 2, 3, 1, -1);
        quat_qf_eval(&result, &qf, &coord);
        print_ibz_hex("qf_eval", &result);

        ibz_vec_4_finalize(&coord);
        ibz_mat_4x4_finalize(&qf);
        ibz_finalize(&result);
    }

    /* ================================================================
     * Section 3: integers, Cornacchia
     * ================================================================ */
    printf("# Section 3: integers\n");

    /* 3a: ibz_cornacchia_prime: n=1, p=5 */
    {
        ibz_t x, y, n, p;
        ibz_init(&x); ibz_init(&y); ibz_init(&n); ibz_init(&p);
        ibz_set(&n, 1); ibz_set(&p, 5);
        int ok = ibz_cornacchia_prime(&x, &y, &n, &p);
        printf("corn_1_5_ok=%d\n", ok);
        if (ok) {
            print_ibz_hex("corn_1_5_x", &x);
            print_ibz_hex("corn_1_5_y", &y);
        }
        ibz_finalize(&x); ibz_finalize(&y); ibz_finalize(&n); ibz_finalize(&p);
    }

    /* 3b: ibz_cornacchia_prime: n=1, p=2 */
    {
        ibz_t x, y, n, p;
        ibz_init(&x); ibz_init(&y); ibz_init(&n); ibz_init(&p);
        ibz_set(&n, 1); ibz_set(&p, 2);
        int ok = ibz_cornacchia_prime(&x, &y, &n, &p);
        printf("corn_1_2_ok=%d\n", ok);
        if (ok) {
            print_ibz_hex("corn_1_2_x", &x);
            print_ibz_hex("corn_1_2_y", &y);
        }
        ibz_finalize(&x); ibz_finalize(&y); ibz_finalize(&n); ibz_finalize(&p);
    }

    /* 3c: ibz_cornacchia_prime: n=1, p=41 */
    {
        ibz_t x, y, n, p;
        ibz_init(&x); ibz_init(&y); ibz_init(&n); ibz_init(&p);
        ibz_set(&n, 1); ibz_set(&p, 41);
        int ok = ibz_cornacchia_prime(&x, &y, &n, &p);
        printf("corn_1_41_ok=%d\n", ok);
        if (ok) {
            print_ibz_hex("corn_1_41_x", &x);
            print_ibz_hex("corn_1_41_y", &y);
        }
        ibz_finalize(&x); ibz_finalize(&y); ibz_finalize(&n); ibz_finalize(&p);
    }

    /* 3d: ibz_cornacchia_prime: n=2, p=3 */
    {
        ibz_t x, y, n, p;
        ibz_init(&x); ibz_init(&y); ibz_init(&n); ibz_init(&p);
        ibz_set(&n, 2); ibz_set(&p, 3);
        int ok = ibz_cornacchia_prime(&x, &y, &n, &p);
        printf("corn_2_3_ok=%d\n", ok);
        if (ok) {
            print_ibz_hex("corn_2_3_x", &x);
            print_ibz_hex("corn_2_3_y", &y);
        }
        ibz_finalize(&x); ibz_finalize(&y); ibz_finalize(&n); ibz_finalize(&p);
    }

    /* 3e: ibz_cornacchia_prime: n=3, p=7 */
    {
        ibz_t x, y, n, p;
        ibz_init(&x); ibz_init(&y); ibz_init(&n); ibz_init(&p);
        ibz_set(&n, 3); ibz_set(&p, 7);
        int ok = ibz_cornacchia_prime(&x, &y, &n, &p);
        printf("corn_3_7_ok=%d\n", ok);
        if (ok) {
            print_ibz_hex("corn_3_7_x", &x);
            print_ibz_hex("corn_3_7_y", &y);
        }
        ibz_finalize(&x); ibz_finalize(&y); ibz_finalize(&n); ibz_finalize(&p);
    }

    /* 3f: ibz_cornacchia_prime: n=1, p=7 (no solution) */
    {
        ibz_t x, y, n, p;
        ibz_init(&x); ibz_init(&y); ibz_init(&n); ibz_init(&p);
        ibz_set(&n, 1); ibz_set(&p, 7);
        int ok = ibz_cornacchia_prime(&x, &y, &n, &p);
        printf("corn_1_7_ok=%d\n", ok);
        ibz_finalize(&x); ibz_finalize(&y); ibz_finalize(&n); ibz_finalize(&p);
    }

    /* 3g: ibz_cornacchia_prime with larger prime: n=1, p=104729 */
    {
        ibz_t x, y, n, p;
        ibz_init(&x); ibz_init(&y); ibz_init(&n); ibz_init(&p);
        ibz_set(&n, 1); ibz_set(&p, 104729);
        int ok = ibz_cornacchia_prime(&x, &y, &n, &p);
        printf("corn_1_104729_ok=%d\n", ok);
        if (ok) {
            print_ibz_hex("corn_1_104729_x", &x);
            print_ibz_hex("corn_1_104729_y", &y);
            /* Verify */
            ibz_t xx, yy, check;
            ibz_init(&xx); ibz_init(&yy); ibz_init(&check);
            ibz_mul(&xx, &x, &x);
            ibz_mul(&yy, &y, &y);
            ibz_add(&check, &xx, &yy);
            printf("corn_1_104729_verify=%d\n", ibz_cmp(&check, &p) == 0 ? 1 : 0);
            ibz_finalize(&xx); ibz_finalize(&yy); ibz_finalize(&check);
        }
        ibz_finalize(&x); ibz_finalize(&y); ibz_finalize(&n); ibz_finalize(&p);
    }

    printf("PASS\n");
    return 0;
}
