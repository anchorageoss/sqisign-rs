/*
 * Cross-validation harness for lll.rs, normeq.rs, lat_ball.rs
 *
 * Exercises LLL reduction, basis reduction, lattice O0 construction,
 * and bound parallelogram computation on hardcoded inputs and prints
 * hex-encoded results to stdout.
 *
 * Build: see build_lll.sh
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <assert.h>

#include "intbig_internal.h"
#include "internal.h"
#include "quaternion.h"
#include "lll_internals.h"

static void print_ibz_hex(const char *label, const ibz_t *x) {
    int sz = ibz_size_in_base(x, 16) + 3;
    char buf[sz];
    memset(buf, 0, sz);
    ibz_convert_to_str(x, buf, 16);
    printf("%s=%s\n", label, buf);
}

static void print_vec4(const char *prefix, const ibz_vec_4_t *v) {
    char label[128];
    for (int i = 0; i < 4; i++) {
        snprintf(label, sizeof(label), "%s_%d", prefix, i);
        print_ibz_hex(label, &(*v)[i]);
    }
}

static void print_mat4x4(const char *prefix, const ibz_mat_4x4_t *m) {
    char label[128];
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            snprintf(label, sizeof(label), "%s_%d%d", prefix, i, j);
            print_ibz_hex(label, &(*m)[i][j]);
        }
    }
}

/* Helper: set the standard O0 order */
static void set_O0(quat_lattice_t *O0) {
    for (int i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++)
            ibz_set(&(O0->basis[i][j]), 0);
    ibz_set(&(O0->denom), 2);
    ibz_set(&(O0->basis[0][0]), 2);
    ibz_set(&(O0->basis[1][1]), 2);
    ibz_set(&(O0->basis[2][2]), 1);
    ibz_set(&(O0->basis[1][2]), 1);
    ibz_set(&(O0->basis[3][3]), 1);
    ibz_set(&(O0->basis[0][3]), 1);
}

/* -----------------------------------------------------------------------
 * Test 1: quat_lattice_lll on O0 with p=103
 * ----------------------------------------------------------------------- */
static void test_lattice_lll_p103(void) {
    printf("=== test_lattice_lll_p103 ===\n");
    fflush(stdout);

    quat_alg_t alg;
    ibz_t p;
    ibz_init(&p);
    ibz_set(&p, 103);
    quat_alg_init_set(&alg, &p);

    quat_lattice_t lat;
    quat_lattice_init(&lat);
    set_O0(&lat);

    ibz_mat_4x4_t red;
    ibz_mat_4x4_init(&red);

    int ret = quat_lattice_lll(&red, &lat, &alg);
    printf("lattice_lll_p103_ret=%d\n", ret);
    print_mat4x4("lattice_lll_p103_red", &red);

    ibz_mat_4x4_finalize(&red);
    quat_lattice_finalize(&lat);
    quat_alg_finalize(&alg);
    ibz_finalize(&p);
}

/* -----------------------------------------------------------------------
 * Test 2: quat_lattice_lll on O0 with p=11
 * ----------------------------------------------------------------------- */
static void test_lattice_lll_p11(void) {
    printf("=== test_lattice_lll_p11 ===\n");
    fflush(stdout);

    quat_alg_t alg;
    ibz_t p;
    ibz_init(&p);
    ibz_set(&p, 11);
    quat_alg_init_set(&alg, &p);

    quat_lattice_t lat;
    quat_lattice_init(&lat);
    set_O0(&lat);

    ibz_mat_4x4_t red;
    ibz_mat_4x4_init(&red);

    int ret = quat_lattice_lll(&red, &lat, &alg);
    printf("lattice_lll_p11_ret=%d\n", ret);
    print_mat4x4("lattice_lll_p11_red", &red);

    ibz_mat_4x4_finalize(&red);
    quat_lattice_finalize(&lat);
    quat_alg_finalize(&alg);
    ibz_finalize(&p);
}

/* -----------------------------------------------------------------------
 * Test 3: quat_lll_verify on reduced basis (p=103)
 * ----------------------------------------------------------------------- */
static void test_lll_verify(void) {
    printf("=== test_lll_verify ===\n");
    fflush(stdout);

    quat_alg_t alg;
    ibz_t p;
    ibz_init(&p);
    ibz_set(&p, 103);
    quat_alg_init_set(&alg, &p);

    quat_lattice_t lat;
    quat_lattice_init(&lat);
    set_O0(&lat);

    ibz_mat_4x4_t red;
    ibz_mat_4x4_init(&red);
    quat_lattice_lll(&red, &lat, &alg);

    /* Set up the delta and eta parameters */
    ibq_t delta, eta;
    ibq_init(&delta);
    ibq_init(&eta);

    ibz_t num, den;
    ibz_init(&num);
    ibz_init(&den);

    /* delta = 99/100 */
    ibz_set(&num, 99);
    ibz_set(&den, 100);
    ibq_set(&delta, &num, &den);

    /* eta = 51/100 */
    ibz_set(&num, 51);
    ibz_set(&den, 100);
    ibq_set(&eta, &num, &den);

    int verified = quat_lll_verify(&red, &delta, &eta, &alg);
    printf("lll_verify_result=%d\n", verified);

    ibz_finalize(&num);
    ibz_finalize(&den);
    ibq_finalize(&delta);
    ibq_finalize(&eta);
    ibz_mat_4x4_finalize(&red);
    quat_lattice_finalize(&lat);
    quat_alg_finalize(&alg);
    ibz_finalize(&p);
}

/* -----------------------------------------------------------------------
 * Test 4: quat_lideal_reduce_basis
 * ----------------------------------------------------------------------- */
static void test_reduce_basis(void) {
    printf("=== test_reduce_basis ===\n");
    fflush(stdout);

    quat_alg_t alg;
    ibz_t p;
    ibz_init(&p);
    ibz_set(&p, 103);
    quat_alg_init_set(&alg, &p);

    /* Create a parent order lattice */
    quat_lattice_t parent;
    quat_lattice_init(&parent);
    set_O0(&parent);

    quat_left_ideal_t lideal;
    quat_left_ideal_init(&lideal);

    /* Set up the ideal's lattice */
    set_O0(&lideal.lattice);
    lideal.parent_order = &parent;

    /* Set norm = 1 */
    ibz_set(&lideal.norm, 1);

    ibz_mat_4x4_t reduced, gram;
    ibz_mat_4x4_init(&reduced);
    ibz_mat_4x4_init(&gram);

    quat_lideal_reduce_basis(&reduced, &gram, &lideal, &alg);

    print_mat4x4("reduce_basis_red", &reduced);
    print_mat4x4("reduce_basis_gram", &gram);

    ibz_mat_4x4_finalize(&reduced);
    ibz_mat_4x4_finalize(&gram);
    quat_left_ideal_finalize(&lideal);
    quat_lattice_finalize(&parent);
    quat_alg_finalize(&alg);
    ibz_finalize(&p);
}

/* -----------------------------------------------------------------------
 * Test 5: quat_lattice_O0_set
 * ----------------------------------------------------------------------- */
static void test_lattice_o0_set(void) {
    printf("=== test_lattice_o0_set ===\n");
    fflush(stdout);

    quat_lattice_t O0;
    quat_lattice_init(&O0);

    quat_lattice_O0_set(&O0);

    print_ibz_hex("o0_denom", &O0.denom);
    print_mat4x4("o0_basis", &O0.basis);

    quat_lattice_finalize(&O0);
}

/* -----------------------------------------------------------------------
 * Test 6: quat_change_to_O0_basis
 * ----------------------------------------------------------------------- */
static void test_change_to_o0_basis(void) {
    printf("=== test_change_to_o0_basis ===\n");
    fflush(stdout);

    quat_alg_elem_t el;
    quat_alg_elem_init(&el);
    ibz_set(&el.denom, 2);
    ibz_set(&el.coord[0], 2);
    ibz_set(&el.coord[1], 7);
    ibz_set(&el.coord[2], 1);
    ibz_set(&el.coord[3], -4);

    ibz_vec_4_t vec;
    ibz_vec_4_init(&vec);

    quat_change_to_O0_basis(&vec, &el);

    print_vec4("o0_basis_vec", &vec);

    ibz_vec_4_finalize(&vec);
    quat_alg_elem_finalize(&el);
}

/* -----------------------------------------------------------------------
 * Test 7: quat_lattice_bound_parallelogram
 * ----------------------------------------------------------------------- */
static void test_bound_parallelogram(void) {
    printf("=== test_bound_parallelogram ===\n");
    fflush(stdout);

    /* Gram matrix: identity * 4 */
    ibz_mat_4x4_t G;
    ibz_mat_4x4_init(&G);
    for (int i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++)
            ibz_set(&G[i][j], i == j ? 4 : 0);

    ibz_t radius;
    ibz_init(&radius);
    ibz_set(&radius, 100);

    ibz_vec_4_t box;
    ibz_vec_4_init(&box);

    ibz_mat_4x4_t U;
    ibz_mat_4x4_init(&U);

    int ret = quat_lattice_bound_parallelogram(&box, &U, &G, &radius);
    printf("bound_para_ret=%d\n", ret);
    print_vec4("bound_para_box", &box);
    print_mat4x4("bound_para_U", &U);

    ibz_mat_4x4_finalize(&U);
    ibz_vec_4_finalize(&box);
    ibz_finalize(&radius);
    ibz_mat_4x4_finalize(&G);
}

int main(void) {
    test_lattice_lll_p103();
    test_lattice_lll_p11();
    test_lll_verify();
    test_reduce_basis();
    test_lattice_o0_set();
    test_change_to_o0_basis();
    test_bound_parallelogram();
    return 0;
}
