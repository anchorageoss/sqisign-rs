/*
 * Cross-validation harness for intbig.c → sqisign-quaternion::intbig.
 *
 * Exercises big integer operations on hardcoded inputs and prints
 * hex-encoded results to stdout. The Rust test
 * (tests/c_crossvalidate_intbig.rs) performs the identical operations
 * and compares byte-for-byte.
 *
 * Build: see build_intbig.sh
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

/* Must come before intbig.h includes, which relies on tutil.h needing RADIX */
#include "intbig_internal.h"

/* Print an ibz_t as hex using ibz_convert_to_str */
static void print_ibz_hex(const char *label, const ibz_t *x) {
    int sz = ibz_size_in_base(x, 16) + 3;
    char buf[sz];
    memset(buf, 0, sz);
    ibz_convert_to_str(x, buf, 16);
    printf("%s=%s\n", label, buf);
}

int main(void) {
    ibz_t a, b, c, q, r, g, inv, s;
    ibz_init(&a);
    ibz_init(&b);
    ibz_init(&c);
    ibz_init(&q);
    ibz_init(&r);
    ibz_init(&g);
    ibz_init(&inv);
    ibz_init(&s);

    /* ------ Section 1: Basic arithmetic ------ */
    printf("# Section 1: Basic arithmetic\n");

    ibz_set_from_str(&a, "deadbeef12345678cafebabe", 16);
    ibz_set_from_str(&b, "1111111122222222333333334444444455555555", 16);

    /* add */
    ibz_add(&c, &a, &b);
    print_ibz_hex("add", &c);

    /* sub */
    ibz_sub(&c, &a, &b);
    print_ibz_hex("sub", &c);

    /* mul */
    ibz_mul(&c, &a, &b);
    print_ibz_hex("mul", &c);

    /* neg */
    ibz_neg(&c, &a);
    print_ibz_hex("neg_a", &c);

    /* abs of negative */
    ibz_abs(&c, &c);
    print_ibz_hex("abs_neg_a", &c);

    /* ------ Section 2: Division ------ */
    printf("# Section 2: Division\n");

    /* Truncated division (ibz_div = mpz_tdiv_qr) */
    ibz_set_from_str(&a, "aaaaaaaabbbbbbbbccccccccdddddddd", 16);
    ibz_set_from_str(&b, "1111111122222222", 16);
    ibz_div(&q, &r, &a, &b);
    print_ibz_hex("tdiv_q", &q);
    print_ibz_hex("tdiv_r", &r);

    /* Truncated division with negative dividend */
    ibz_neg(&a, &a);
    ibz_div(&q, &r, &a, &b);
    print_ibz_hex("tdiv_negdiv_q", &q);
    print_ibz_hex("tdiv_negdiv_r", &r);
    ibz_neg(&a, &a); /* restore */

    /* Floor division */
    ibz_neg(&a, &a);
    ibz_div_floor(&q, &r, &a, &b);
    print_ibz_hex("fdiv_negdiv_q", &q);
    print_ibz_hex("fdiv_negdiv_r", &r);
    ibz_neg(&a, &a); /* restore */

    /* ibz_mod: always non-negative */
    ibz_set_from_str(&a, "-deadbeefcafebabe", 16);
    ibz_set_from_str(&b, "1234567890abcdef", 16);
    ibz_mod(&r, &a, &b);
    print_ibz_hex("mod_neg", &r);

    /* div_2exp */
    ibz_set_from_str(&a, "ffffffffffffffffffffffffffffffff", 16);
    ibz_div_2exp(&q, &a, 64);
    print_ibz_hex("div2exp_64", &q);

    /* div_2exp negative */
    ibz_neg(&a, &a);
    ibz_div_2exp(&q, &a, 17);
    print_ibz_hex("div2exp_neg_17", &q);
    ibz_neg(&a, &a); /* restore */

    /* ------ Section 3: Number theory ------ */
    printf("# Section 3: Number theory\n");

    /* GCD */
    ibz_set_from_str(&a, "3b9aca00", 16);  /* 1000000000 */
    ibz_set_from_str(&b, "e8d4a51000", 16); /* 1000000000000 */
    ibz_gcd(&g, &a, &b);
    print_ibz_hex("gcd", &g);

    /* pow */
    ibz_set_from_str(&a, "ff", 16);
    ibz_pow(&c, &a, 7);
    print_ibz_hex("pow_ff_7", &c);

    /* pow_mod */
    ibz_set_from_str(&a, "deadbeef", 16);
    ibz_set_from_str(&b, "1234567890", 16);
    ibz_set_from_str(&c, "ffffffffffffffc5", 16); /* = 2^64 - 59, a prime */
    ibz_pow_mod(&q, &a, &b, &c);
    print_ibz_hex("powmod", &q);

    /* invmod */
    ibz_set_from_str(&a, "deadbeef12345678", 16);
    ibz_set_from_str(&b, "ffffffffffffffc5", 16);
    int ok = ibz_invmod(&inv, &a, &b);
    printf("invmod_ok=%d\n", ok);
    print_ibz_hex("invmod", &inv);

    /* Legendre */
    ibz_set_from_str(&a, "3", 16);
    ibz_set_from_str(&b, "ffffffffffffffc5", 16); /* prime */
    int leg = ibz_legendre(&a, &b);
    printf("legendre_3=%d\n", leg);

    ibz_set(&a, 2);
    leg = ibz_legendre(&a, &b);
    printf("legendre_2=%d\n", leg);

    /* two_adic */
    ibz_set_from_str(&a, "abcdef0000000000", 16);
    int tav = ibz_two_adic(&a);
    printf("two_adic=%d\n", tav);

    /* sqrt_floor */
    ibz_set_from_str(&a, "10000000000000000", 16); /* 2^64 */
    ibz_sqrt_floor(&s, &a);
    print_ibz_hex("sqrt_floor", &s);

    /* sqrt (perfect square) */
    ibz_set_from_str(&a, "1000000000000", 16); /* 2^48 */
    ibz_mul(&a, &a, &a); /* 2^96 */
    ok = ibz_sqrt(&s, &a);
    printf("sqrt_perfect_ok=%d\n", ok);
    print_ibz_hex("sqrt_perfect", &s);

    /* sqrt_mod_p: p = 2^64 - 59 (prime, p ≡ 1 mod 8 → Tonelli-Shanks) */
    ibz_set_from_str(&b, "ffffffffffffffc5", 16);
    ibz_set_from_str(&a, "9", 16);
    ok = ibz_sqrt_mod_p(&s, &a, &b);
    printf("sqrtmodp_9_ok=%d\n", ok);
    /* Verify: s^2 mod p == 9 */
    ibz_pow_mod(&c, &s, &ibz_const_two, &b);
    print_ibz_hex("sqrtmodp_9_sq", &c);

    /* sqrt_mod_p: p ≡ 3 mod 4 */
    ibz_set_from_str(&b, "1f", 16); /* 31, prime, 31 % 4 == 3 */
    ibz_set_from_str(&a, "4", 16);
    ok = ibz_sqrt_mod_p(&s, &a, &b);
    printf("sqrtmodp_p3m4_ok=%d\n", ok);
    ibz_pow_mod(&c, &s, &ibz_const_two, &b);
    print_ibz_hex("sqrtmodp_p3m4_sq", &c);

    /* sqrt_mod_p: p ≡ 5 mod 8 */
    ibz_set_from_str(&b, "d", 16); /* 13, prime, 13 % 8 == 5 */
    ibz_set_from_str(&a, "4", 16);
    ok = ibz_sqrt_mod_p(&s, &a, &b);
    printf("sqrtmodp_p5m8_ok=%d\n", ok);
    ibz_pow_mod(&c, &s, &ibz_const_two, &b);
    print_ibz_hex("sqrtmodp_p5m8_sq", &c);

    /* ------ Section 4: Digit conversion ------ */
    printf("# Section 4: Digit conversion\n");

    /* copy_digits: import 2 limbs [0x0002, 0x0001] = 1 * 2^64 + 2 */
    digit_t limbs_in[2] = { 0x0000000000000002ULL, 0x0000000000000001ULL };
    ibz_copy_digits(&a, limbs_in, 2);
    print_ibz_hex("copy_digits_2_1", &a);

    /* to_digits: export and re-import */
    digit_t limbs_out[2] = { 0, 0 };
    ibz_to_digits(limbs_out, &a);
    printf("to_digits_0=%016llx\n", (unsigned long long)limbs_out[0]);
    printf("to_digits_1=%016llx\n", (unsigned long long)limbs_out[1]);

    /* ------ Section 5: Comparison / predicates ------ */
    printf("# Section 5: Comparison\n");

    ibz_set(&a, 0);
    printf("is_zero_0=%d\n", ibz_is_zero(&a));
    ibz_set(&a, 1);
    printf("is_one_1=%d\n", ibz_is_one(&a));
    printf("is_even_1=%d\n", ibz_is_even(&a));
    printf("is_odd_1=%d\n", ibz_is_odd(&a));
    ibz_set(&a, 42);
    printf("is_even_42=%d\n", ibz_is_even(&a));
    printf("is_odd_42=%d\n", ibz_is_odd(&a));
    printf("bitsize_42=%d\n", ibz_bitsize(&a));

    ibz_set_from_str(&a, "deadbeef12345678", 16);
    printf("get_lo32=0x%x\n", (unsigned)ibz_get(&a));

    ibz_set(&a, 7);
    ibz_set(&b, 7);
    printf("cmp_eq=%d\n", ibz_cmp(&a, &b) == 0 ? 1 : 0);
    ibz_set(&a, 3);
    printf("cmp_lt=%d\n", ibz_cmp(&a, &b) < 0 ? 1 : 0);
    ibz_set(&a, 10);
    printf("cmp_gt=%d\n", ibz_cmp(&a, &b) > 0 ? 1 : 0);

    /* mod_ui */
    ibz_set_from_str(&a, "2113309833171849999003363", 10);
    printf("mod_ui_3=%lu\n", ibz_mod_ui(&a, 3));
    printf("mod_ui_2=%lu\n", ibz_mod_ui(&a, 2));

    /* divides */
    ibz_set(&a, 12);
    ibz_set(&b, 3);
    printf("divides_12_3=%d\n", ibz_divides(&a, &b));
    ibz_set(&b, 5);
    printf("divides_12_5=%d\n", ibz_divides(&a, &b));

    /* probab_prime */
    ibz_set(&a, 17);
    printf("prime_17=%d\n", ibz_probab_prime(&a, 25) > 0 ? 1 : 0);
    ibz_set(&a, 15);
    printf("prime_15=%d\n", ibz_probab_prime(&a, 25) > 0 ? 1 : 0);

    ibz_finalize(&a);
    ibz_finalize(&b);
    ibz_finalize(&c);
    ibz_finalize(&q);
    ibz_finalize(&r);
    ibz_finalize(&g);
    ibz_finalize(&inv);
    ibz_finalize(&s);

    printf("PASS\n");
    return 0;
}
