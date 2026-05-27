/*
 * Cross-validation harness for algebra.c, lattice.c, ideal.c, hnf.c
 * → sqisign-quaternion::{algebra, lattice, ideal, hnf}.
 *
 * Exercises quaternion algebra, lattice, and ideal operations on hardcoded
 * inputs and prints hex-encoded results to stdout. The Rust test
 * (tests/c_crossvalidate_quat.rs) performs the identical operations
 * and compares byte-for-byte.
 *
 * Build: see build_quat.sh
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <assert.h>

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

static void print_elem(const char *prefix, const quat_alg_elem_t *e) {
    char label[128];
    snprintf(label, sizeof(label), "%s_denom", prefix);
    print_ibz_hex(label, &e->denom);
    for (int i = 0; i < 4; i++) {
        snprintf(label, sizeof(label), "%s_coord%d", prefix, i);
        print_ibz_hex(label, &e->coord[i]);
    }
}

static void print_lattice(const char *prefix, const quat_lattice_t *lat) {
    char label[128];
    snprintf(label, sizeof(label), "%s_denom", prefix);
    print_ibz_hex(label, &lat->denom);
    snprintf(label, sizeof(label), "%s_basis", prefix);
    print_mat4x4(label, &lat->basis);
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

int main(void) {
    /* ================================================================
     * Section 1: algebra, coord_mul, add, sub, mul, norm, conj, normalize
     * ================================================================ */
    printf("# Section 1: algebra\n");

    /* 1a: coord_mul */
    {
        quat_alg_t alg;
        quat_alg_init_set_ui(&alg, 7);
        ibz_vec_4_t a, b, c;
        ibz_vec_4_init(&a);
        ibz_vec_4_init(&b);
        ibz_vec_4_init(&c);
        ibz_set(&a[0], 152); ibz_set(&a[1], 57);
        ibz_set(&a[2], 190); ibz_set(&a[3], 28);
        ibz_set(&b[0], 165); ibz_set(&b[1], 35);
        ibz_set(&b[2], 231); ibz_set(&b[3], 770);
        quat_alg_coord_mul(&c, &a, &b, &alg);
        print_vec4("coord_mul_p7", &c);

        ibz_set(&alg.p, 11);
        quat_alg_coord_mul(&c, &a, &b, &alg);
        print_vec4("coord_mul_p11", &c);
        ibz_vec_4_finalize(&a);
        ibz_vec_4_finalize(&b);
        ibz_vec_4_finalize(&c);
        quat_alg_finalize(&alg);
    }

    /* 1b: mul with denom */
    {
        quat_alg_t alg;
        quat_alg_init_set_ui(&alg, 7);
        quat_alg_elem_t a, b, c;
        quat_alg_elem_init(&a);
        quat_alg_elem_init(&b);
        quat_alg_elem_init(&c);
        ibz_set(&a.coord[0], 152); ibz_set(&a.coord[1], 57);
        ibz_set(&a.coord[2], 190); ibz_set(&a.coord[3], 28);
        ibz_set(&a.denom, 76);
        ibz_set(&b.coord[0], 165); ibz_set(&b.coord[1], 35);
        ibz_set(&b.coord[2], 231); ibz_set(&b.coord[3], 770);
        ibz_set(&b.denom, 385);
        quat_alg_mul(&c, &a, &b, &alg);
        print_elem("mul_p7", &c);
        quat_alg_elem_finalize(&a);
        quat_alg_elem_finalize(&b);
        quat_alg_elem_finalize(&c);
        quat_alg_finalize(&alg);
    }

    /* 1c: add */
    {
        quat_alg_elem_t a, b, c;
        quat_alg_elem_init(&a);
        quat_alg_elem_init(&b);
        quat_alg_elem_init(&c);
        ibz_set(&a.coord[0], -12); ibz_set(&a.coord[1], 0);
        ibz_set(&a.coord[2], -7); ibz_set(&a.coord[3], 19);
        ibz_set(&a.denom, 9);
        ibz_set(&b.coord[0], -6); ibz_set(&b.coord[1], 2);
        ibz_set(&b.coord[2], 7); ibz_set(&b.coord[3], -19);
        ibz_set(&b.denom, 3);
        quat_alg_add(&c, &a, &b);
        print_elem("add_1", &c);
        ibz_set(&b.denom, 6);
        quat_alg_add(&c, &a, &b);
        print_elem("add_2", &c);
        quat_alg_elem_finalize(&a);
        quat_alg_elem_finalize(&b);
        quat_alg_elem_finalize(&c);
    }

    /* 1d: sub */
    {
        quat_alg_elem_t a, b, c;
        quat_alg_elem_init(&a);
        quat_alg_elem_init(&b);
        quat_alg_elem_init(&c);
        ibz_set(&a.coord[0], -12); ibz_set(&a.coord[1], 0);
        ibz_set(&a.coord[2], -7); ibz_set(&a.coord[3], 19);
        ibz_set(&a.denom, 9);
        ibz_set(&b.coord[0], -6); ibz_set(&b.coord[1], 2);
        ibz_set(&b.coord[2], 7); ibz_set(&b.coord[3], -19);
        ibz_set(&b.denom, 3);
        quat_alg_sub(&c, &a, &b);
        print_elem("sub_1", &c);
        quat_alg_elem_finalize(&a);
        quat_alg_elem_finalize(&b);
        quat_alg_elem_finalize(&c);
    }

    /* 1e: norm */
    {
        quat_alg_t alg;
        quat_alg_init_set_ui(&alg, 11);
        quat_alg_elem_t a;
        ibz_t num, denom;
        quat_alg_elem_init(&a);
        ibz_init(&num); ibz_init(&denom);

        ibz_set(&a.coord[0], 1); ibz_set(&a.coord[1], 5);
        ibz_set(&a.coord[2], 7); ibz_set(&a.coord[3], 2);
        ibz_set(&a.denom, 2);
        quat_alg_norm(&num, &denom, &a, &alg);
        print_ibz_hex("norm1_num", &num);
        print_ibz_hex("norm1_denom", &denom);

        ibz_set(&a.coord[0], 152); ibz_set(&a.coord[1], 57);
        ibz_set(&a.coord[2], 190); ibz_set(&a.coord[3], 28);
        ibz_set(&a.denom, 76);
        quat_alg_norm(&num, &denom, &a, &alg);
        print_ibz_hex("norm2_num", &num);
        print_ibz_hex("norm2_denom", &denom);

        ibz_set(&alg.p, 7);
        quat_alg_norm(&num, &denom, &a, &alg);
        print_ibz_hex("norm3_num", &num);
        print_ibz_hex("norm3_denom", &denom);

        quat_alg_elem_finalize(&a);
        ibz_finalize(&num); ibz_finalize(&denom);
        quat_alg_finalize(&alg);
    }

    /* 1f: conj */
    {
        quat_alg_elem_t a, c;
        quat_alg_elem_init(&a);
        quat_alg_elem_init(&c);
        ibz_set(&a.coord[0], -125); ibz_set(&a.coord[1], 2);
        ibz_set(&a.coord[2], 0); ibz_set(&a.coord[3], -30);
        ibz_set(&a.denom, 25);
        quat_alg_conj(&c, &a);
        print_elem("conj_1", &c);
        quat_alg_elem_finalize(&a);
        quat_alg_elem_finalize(&c);
    }

    /* 1g: normalize */
    {
        quat_alg_elem_t x;
        quat_alg_elem_init(&x);
        ibz_set(&x.coord[0], -36); ibz_set(&x.coord[1], 18);
        ibz_set(&x.coord[2], 0); ibz_set(&x.coord[3], -300);
        ibz_set(&x.denom, 48);
        quat_alg_normalize(&x);
        print_elem("normalize_1", &x);

        ibz_set(&x.coord[0], -36); ibz_set(&x.coord[1], 18);
        ibz_set(&x.coord[2], 0); ibz_set(&x.coord[3], -300);
        ibz_set(&x.denom, -6);
        quat_alg_normalize(&x);
        print_elem("normalize_2", &x);
        quat_alg_elem_finalize(&x);
    }

    /* ================================================================
     * Section 2: HNF
     * ================================================================ */
    printf("# Section 2: hnf\n");

    /* 2a: hnf of a non-trivial lattice */
    {
        quat_lattice_t lat;
        quat_lattice_init(&lat);
        for (int i = 0; i < 4; i++)
            for (int j = 0; j < 4; j++)
                ibz_set(&lat.basis[i][j], 0);
        ibz_set(&lat.basis[0][0], 1);
        ibz_set(&lat.basis[0][3], -1);
        ibz_set(&lat.basis[1][1], -2);
        ibz_set(&lat.basis[2][2], 1);
        ibz_set(&lat.basis[2][1], 1);
        ibz_set(&lat.basis[3][3], -3);
        ibz_set(&lat.denom, 6);
        quat_lattice_hnf(&lat);
        print_lattice("hnf_1", &lat);
        quat_lattice_finalize(&lat);
    }

    /* ================================================================
     * Section 3: lattice operations
     * ================================================================ */
    printf("# Section 3: lattice\n");

    /* 3a: lattice_add */
    {
        quat_lattice_t lat1, lat2, sum;
        quat_lattice_init(&lat1);
        quat_lattice_init(&lat2);
        quat_lattice_init(&sum);
        ibz_mat_4x4_zero(&lat1.basis);
        ibz_mat_4x4_zero(&lat2.basis);
        ibz_set(&lat1.basis[0][0], 4);
        ibz_set(&lat1.basis[0][2], 3);
        ibz_set(&lat2.basis[0][0], 1);
        ibz_set(&lat2.basis[0][3], -1);
        ibz_set(&lat1.basis[1][1], 5);
        ibz_set(&lat2.basis[1][1], -2);
        ibz_set(&lat1.basis[2][2], 3);
        ibz_set(&lat2.basis[2][2], 1);
        ibz_set(&lat2.basis[2][1], 1);
        ibz_set(&lat1.basis[3][3], 7);
        ibz_set(&lat2.basis[3][3], -3);
        ibz_set(&lat1.denom, 4);
        ibz_set(&lat2.denom, 6);
        quat_lattice_add(&sum, &lat1, &lat2);
        print_lattice("lattice_add", &sum);
        quat_lattice_finalize(&lat1);
        quat_lattice_finalize(&lat2);
        quat_lattice_finalize(&sum);
    }

    /* 3b: lattice_intersect */
    {
        quat_lattice_t lat1, lat2, inter;
        quat_lattice_init(&lat1);
        quat_lattice_init(&lat2);
        quat_lattice_init(&inter);
        ibz_mat_4x4_zero(&lat1.basis);
        ibz_mat_4x4_zero(&lat2.basis);
        ibz_set(&lat1.basis[0][0], 4);
        ibz_set(&lat1.basis[0][2], 3);
        ibz_set(&lat2.basis[0][0], 1);
        ibz_set(&lat2.basis[0][3], -1);
        ibz_set(&lat1.basis[1][1], 5);
        ibz_set(&lat2.basis[1][1], -2);
        ibz_set(&lat1.basis[2][2], 3);
        ibz_set(&lat2.basis[2][2], 1);
        ibz_set(&lat2.basis[2][1], 1);
        ibz_set(&lat1.basis[3][3], 7);
        ibz_set(&lat2.basis[3][3], -3);
        ibz_set(&lat1.denom, 4);
        ibz_set(&lat2.denom, 6);
        quat_lattice_hnf(&lat1);
        quat_lattice_hnf(&lat2);
        quat_lattice_intersect(&inter, &lat1, &lat2);
        print_lattice("lattice_inter", &inter);
        quat_lattice_finalize(&lat1);
        quat_lattice_finalize(&lat2);
        quat_lattice_finalize(&inter);
    }

    /* 3c: lattice_index */
    {
        quat_lattice_t sublat, overlat;
        ibz_t index;
        ibz_init(&index);
        quat_lattice_init(&sublat);
        quat_lattice_init(&overlat);
        ibz_mat_4x4_zero(&sublat.basis);
        ibz_mat_4x4_identity(&overlat.basis);
        ibz_set(&overlat.denom, 2);
        ibz_set(&sublat.basis[0][0], 2);
        ibz_set(&sublat.basis[0][2], 1);
        ibz_set(&sublat.basis[1][1], 4);
        ibz_set(&sublat.basis[1][2], 2);
        ibz_set(&sublat.basis[1][3], 3);
        ibz_set(&sublat.basis[2][2], 1);
        ibz_set(&sublat.basis[3][3], 1);
        ibz_set(&sublat.denom, 2);
        quat_lattice_index(&index, &sublat, &overlat);
        print_ibz_hex("lattice_index", &index);
        ibz_finalize(&index);
        quat_lattice_finalize(&sublat);
        quat_lattice_finalize(&overlat);
    }

    /* 3d: lattice_contains */
    {
        quat_lattice_t lat;
        quat_alg_elem_t x;
        ibz_vec_4_t coord;
        quat_lattice_init(&lat);
        quat_alg_elem_init(&x);
        ibz_vec_4_init(&coord);
        ibz_mat_4x4_zero(&lat.basis);
        ibz_set(&lat.basis[0][0], 1);
        ibz_set(&lat.basis[0][3], -1);
        ibz_set(&lat.basis[1][1], -2);
        ibz_set(&lat.basis[2][2], 1);
        ibz_set(&lat.basis[2][1], 1);
        ibz_set(&lat.basis[3][3], -3);
        ibz_set(&lat.denom, 6);
        quat_lattice_hnf(&lat);
        ibz_set(&x.denom, 3);
        ibz_set(&x.coord[0], 1);
        ibz_set(&x.coord[1], -2);
        ibz_set(&x.coord[2], 26);
        ibz_set(&x.coord[3], 9);
        int ok = quat_lattice_contains(&coord, &lat, &x);
        printf("lattice_contains_ok=%d\n", ok);
        print_vec4("lattice_contains_coord", &coord);
        quat_lattice_finalize(&lat);
        quat_alg_elem_finalize(&x);
        ibz_vec_4_finalize(&coord);
    }

    /* 3e: lattice_conjugate_without_hnf */
    {
        quat_lattice_t lat, conj;
        quat_lattice_init(&lat);
        quat_lattice_init(&conj);
        ibz_mat_4x4_zero(&lat.basis);
        ibz_set(&lat.basis[0][0], 4);
        ibz_set(&lat.basis[0][3], 1);
        ibz_set(&lat.basis[1][1], -2);
        ibz_set(&lat.basis[2][2], -1);
        ibz_set(&lat.basis[2][1], -1);
        ibz_set(&lat.basis[3][3], -3);
        ibz_set(&lat.denom, 6);
        quat_lattice_hnf(&lat);
        quat_lattice_conjugate_without_hnf(&conj, &lat);
        quat_lattice_hnf(&conj);
        print_lattice("lattice_conj", &conj);
        quat_lattice_finalize(&lat);
        quat_lattice_finalize(&conj);
    }

    /* 3f: lattice_dual_without_hnf */
    {
        quat_lattice_t lat, dual;
        quat_lattice_init(&lat);
        quat_lattice_init(&dual);
        ibz_mat_4x4_zero(&lat.basis);
        ibz_set(&lat.basis[0][0], 1);
        ibz_set(&lat.basis[0][3], -1);
        ibz_set(&lat.basis[1][1], -2);
        ibz_set(&lat.basis[2][2], 1);
        ibz_set(&lat.basis[2][1], 1);
        ibz_set(&lat.basis[3][3], -3);
        ibz_set(&lat.denom, 6);
        quat_lattice_hnf(&lat);
        quat_lattice_dual_without_hnf(&dual, &lat);
        quat_lattice_hnf(&dual);
        print_lattice("lattice_dual", &dual);
        quat_lattice_finalize(&lat);
        quat_lattice_finalize(&dual);
    }

    /* 3g: lattice_mul */
    {
        quat_lattice_t lat1, lat2, prod;
        quat_alg_t alg;
        quat_lattice_init(&lat1);
        quat_lattice_init(&lat2);
        quat_lattice_init(&prod);
        quat_alg_init_set_ui(&alg, 19);
        ibz_mat_4x4_zero(&lat1.basis);
        ibz_mat_4x4_zero(&lat2.basis);
        ibz_set(&lat1.basis[0][0], 44);
        ibz_set(&lat1.basis[0][2], 3);
        ibz_set(&lat1.basis[0][3], 32);
        ibz_set(&lat2.basis[0][0], 1);
        ibz_set(&lat1.basis[1][1], 5);
        ibz_set(&lat2.basis[1][1], 2);
        ibz_set(&lat1.basis[2][2], 3);
        ibz_set(&lat2.basis[2][2], 1);
        ibz_set(&lat1.basis[3][3], 1);
        ibz_set(&lat2.basis[3][3], 3);
        ibz_set(&lat1.denom, 4);
        ibz_set(&lat2.denom, 6);
        quat_lattice_mul(&prod, &lat1, &lat2, &alg);
        print_lattice("lattice_mul", &prod);
        quat_lattice_finalize(&lat1);
        quat_lattice_finalize(&lat2);
        quat_lattice_finalize(&prod);
        quat_alg_finalize(&alg);
    }

    /* 3h: lattice_gram */
    {
        quat_lattice_t lat;
        ibz_mat_4x4_t gram;
        quat_alg_t alg;
        quat_lattice_init(&lat);
        ibz_mat_4x4_init(&gram);
        quat_alg_init_set_ui(&alg, 103);
        set_O0(&lat);
        quat_lattice_gram(&gram, &lat, &alg);
        print_mat4x4("lattice_gram", &gram);
        quat_lattice_finalize(&lat);
        ibz_mat_4x4_finalize(&gram);
        quat_alg_finalize(&alg);
    }

    /* ================================================================
     * Section 4: ideal operations
     * ================================================================ */
    printf("# Section 4: ideal\n");

    /* 4a: lideal_create_principal with p=367 */
    {
        quat_alg_t alg;
        quat_lattice_t order;
        quat_alg_elem_t gamma;
        quat_left_ideal_t I;
        quat_alg_init_set_ui(&alg, 367);
        quat_lattice_init(&order);
        quat_alg_elem_init(&gamma);
        quat_left_ideal_init(&I);
        set_O0(&order);
        ibz_set(&gamma.coord[0], 219);
        ibz_set(&gamma.coord[1], 200);
        ibz_set(&gamma.coord[2], 78);
        ibz_set(&gamma.coord[3], -1);
        quat_lideal_create_principal(&I, &gamma, &order, &alg);
        print_ibz_hex("principal_norm", &I.norm);
        print_lattice("principal_lat", &I.lattice);
        quat_alg_finalize(&alg);
        quat_lattice_finalize(&order);
        quat_alg_elem_finalize(&gamma);
        quat_left_ideal_finalize(&I);
    }

    /* 4b: lideal_create with p=367, N=31 */
    {
        quat_alg_t alg;
        quat_lattice_t order;
        quat_alg_elem_t gamma;
        ibz_t N;
        quat_left_ideal_t I;
        quat_alg_init_set_ui(&alg, 367);
        quat_lattice_init(&order);
        quat_alg_elem_init(&gamma);
        ibz_init(&N);
        quat_left_ideal_init(&I);
        set_O0(&order);
        ibz_set(&gamma.coord[0], 219);
        ibz_set(&gamma.coord[1], 200);
        ibz_set(&gamma.coord[2], 78);
        ibz_set(&gamma.coord[3], -1);
        ibz_set(&N, 31);
        quat_lideal_create(&I, &gamma, &N, &order, &alg);
        print_ibz_hex("create_norm", &I.norm);
        print_lattice("create_lat", &I.lattice);
        quat_alg_finalize(&alg);
        quat_lattice_finalize(&order);
        quat_alg_elem_finalize(&gamma);
        ibz_finalize(&N);
        quat_left_ideal_finalize(&I);
    }

    /* 4c: lideal_add, lideal_inter, lideal_equals with p=103 */
    {
        quat_alg_t alg;
        quat_lattice_t order;
        quat_alg_elem_t gen1, gen2;
        ibz_t N1, N2;
        quat_left_ideal_t lideal1, lideal2, sum, inter;
        quat_alg_init_set_ui(&alg, 103);
        quat_lattice_init(&order);
        quat_alg_elem_init(&gen1);
        quat_alg_elem_init(&gen2);
        ibz_init(&N1); ibz_init(&N2);
        quat_left_ideal_init(&lideal1);
        quat_left_ideal_init(&lideal2);
        quat_left_ideal_init(&sum);
        quat_left_ideal_init(&inter);
        set_O0(&order);

        ibz_set(&gen1.coord[0], 3); ibz_set(&gen1.coord[1], 5);
        ibz_set(&gen1.coord[2], 7); ibz_set(&gen1.coord[3], 11);
        ibz_set(&N1, 17);
        quat_lideal_create(&lideal1, &gen1, &N1, &order, &alg);

        ibz_set(&gen2.coord[0], -2); ibz_set(&gen2.coord[1], 13);
        ibz_set(&gen2.coord[2], -17); ibz_set(&gen2.coord[3], 19);
        ibz_set(&N2, 43);
        quat_lideal_create(&lideal2, &gen2, &N2, &order, &alg);

        /* sum should be the whole order */
        quat_lideal_add(&sum, &lideal1, &lideal2, &alg);
        print_ibz_hex("add_norm", &sum.norm);
        print_lattice("add_lat", &sum.lattice);

        /* self-intersection */
        quat_lideal_inter(&inter, &lideal1, &lideal1, &alg);
        int eq = quat_lideal_equals(&inter, &lideal1, &alg);
        printf("selfinter_eq=%d\n", eq);
        print_ibz_hex("selfinter_norm", &inter.norm);

        quat_alg_finalize(&alg);
        quat_lattice_finalize(&order);
        quat_alg_elem_finalize(&gen1);
        quat_alg_elem_finalize(&gen2);
        ibz_finalize(&N1); ibz_finalize(&N2);
        quat_left_ideal_finalize(&lideal1);
        quat_left_ideal_finalize(&lideal2);
        quat_left_ideal_finalize(&sum);
        quat_left_ideal_finalize(&inter);
    }

    /* 4d: order_discriminant and order_is_maximal with p=43 */
    {
        quat_alg_t alg;
        quat_lattice_t order;
        ibz_t disc;
        quat_alg_init_set_ui(&alg, 43);
        quat_lattice_init(&order);
        ibz_init(&disc);
        set_O0(&order);
        quat_order_discriminant(&disc, &order, &alg);
        print_ibz_hex("disc_O0", &disc);
        int maximal = quat_order_is_maximal(&order, &alg);
        printf("is_maximal_O0=%d\n", maximal);

        /* Z^4 is not maximal */
        ibz_mat_4x4_identity(&order.basis);
        ibz_set(&order.denom, 1);
        int not_maximal = quat_order_is_maximal(&order, &alg);
        printf("is_maximal_Z4=%d\n", not_maximal);

        ibz_finalize(&disc);
        quat_lattice_finalize(&order);
        quat_alg_finalize(&alg);
    }

    /* 4e: lideal_right_order with p=19 */
    {
        quat_alg_t alg;
        quat_lattice_t order, rorder;
        quat_alg_elem_t gen, test_elem;
        ibz_t norm;
        quat_left_ideal_t lideal;
        quat_alg_init_set_ui(&alg, 19);
        quat_lattice_init(&order);
        quat_lattice_init(&rorder);
        quat_alg_elem_init(&gen);
        quat_alg_elem_init(&test_elem);
        ibz_init(&norm);
        quat_left_ideal_init(&lideal);

        ibz_set(&order.basis[0][0], 4);
        ibz_set(&order.basis[0][1], 0);
        ibz_set(&order.basis[0][2], 2);
        ibz_set(&order.basis[0][3], 2);
        ibz_set(&order.basis[1][0], 0);
        ibz_set(&order.basis[1][1], 8);
        ibz_set(&order.basis[1][2], 4);
        ibz_set(&order.basis[1][3], 3);
        ibz_set(&order.basis[2][0], 0);
        ibz_set(&order.basis[2][1], 0);
        ibz_set(&order.basis[2][2], 2);
        ibz_set(&order.basis[2][3], 0);
        ibz_set(&order.basis[3][0], 0);
        ibz_set(&order.basis[3][1], 0);
        ibz_set(&order.basis[3][2], 0);
        ibz_set(&order.basis[3][3], 1);
        ibz_set(&order.denom, 4);
        quat_alg_elem_set(&gen, 1, 3, 3, 0, 1);
        ibz_set(&norm, 15);
        quat_lideal_create(&lideal, &gen, &norm, &order, &alg);
        quat_lideal_right_order(&rorder, &lideal, &alg);
        print_lattice("rorder", &rorder);

        /* verify it's in HNF */
        int is_hnf = ibz_mat_4x4_is_hnf(&rorder.basis);
        printf("rorder_is_hnf=%d\n", is_hnf);

        /* verify contains 1 */
        quat_alg_elem_set(&test_elem, 1, 1, 0, 0, 0);
        int contains_one = quat_lattice_contains(NULL, &rorder, &test_elem);
        printf("rorder_contains_one=%d\n", contains_one);

        quat_alg_finalize(&alg);
        quat_lattice_finalize(&order);
        quat_lattice_finalize(&rorder);
        quat_alg_elem_finalize(&gen);
        quat_alg_elem_finalize(&test_elem);
        ibz_finalize(&norm);
        quat_left_ideal_finalize(&lideal);
    }

    /* 4f: lideal_inverse_lattice_without_hnf with p=19 */
    {
        quat_alg_t alg;
        quat_lattice_t order, inv, prod;
        quat_alg_elem_t gen;
        ibz_t norm;
        quat_left_ideal_t lideal;
        quat_alg_init_set_ui(&alg, 19);
        quat_lattice_init(&order);
        quat_lattice_init(&inv);
        quat_lattice_init(&prod);
        quat_alg_elem_init(&gen);
        ibz_init(&norm);
        quat_left_ideal_init(&lideal);

        ibz_set(&order.basis[0][0], 4);
        ibz_set(&order.basis[0][1], 0);
        ibz_set(&order.basis[0][2], 2);
        ibz_set(&order.basis[0][3], 2);
        ibz_set(&order.basis[1][0], 0);
        ibz_set(&order.basis[1][1], 8);
        ibz_set(&order.basis[1][2], 4);
        ibz_set(&order.basis[1][3], 3);
        ibz_set(&order.basis[2][0], 0);
        ibz_set(&order.basis[2][1], 0);
        ibz_set(&order.basis[2][2], 2);
        ibz_set(&order.basis[2][3], 0);
        ibz_set(&order.basis[3][0], 0);
        ibz_set(&order.basis[3][1], 0);
        ibz_set(&order.basis[3][2], 0);
        ibz_set(&order.basis[3][3], 1);
        ibz_set(&order.denom, 4);
        quat_alg_elem_set(&gen, 1, 2, 3, 0, 1);
        ibz_set(&norm, 15);
        quat_lideal_create(&lideal, &gen, &norm, &order, &alg);
        quat_lideal_inverse_lattice_without_hnf(&inv, &lideal, &alg);
        quat_lattice_mul(&prod, &lideal.lattice, &inv, &alg);
        int eq = quat_lattice_equal(&prod, &order);
        printf("inv_prod_eq_order=%d\n", eq);

        quat_alg_finalize(&alg);
        quat_lattice_finalize(&order);
        quat_lattice_finalize(&inv);
        quat_lattice_finalize(&prod);
        quat_alg_elem_finalize(&gen);
        ibz_finalize(&norm);
        quat_left_ideal_finalize(&lideal);
    }

    printf("# Done\n");
    return 0;
}
