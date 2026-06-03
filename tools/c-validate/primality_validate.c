/*
 * primality_validate.c, Cross-validation harness for ibz_probab_prime.
 *
 * Links against the C SQIsign reference (which uses mini-gmp internally)
 * and prints the result of mpz_probab_prime_p for a set of test values.
 * The Rust regression tests must produce identical results.
 *
 * Build:
 *   gcc -o primality_cval primality_validate.c -I/path/to/sqisign/include \
 *       -L/path/to/sqisign/build -lsqisign -lgmp -lm
 *
 * Or with mini-gmp (no system GMP):
 *   gcc -o primality_cval primality_validate.c mini-gmp.c -I. -DUSE_MINI_GMP
 *
 * Usage:
 *   ./primality_cval
 *
 * Output: one line per test value, format:
 *   <decimal_value> <result>
 * where result is 0 (composite), 1 (probably prime), or 2 (certainly prime).
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef USE_MINI_GMP
#include "mini-gmp.h"
#else
#include <gmp.h>
#endif

#define REPS 32

static void test_value(const char *decimal_str) {
    mpz_t n;
    mpz_init(n);
    mpz_set_str(n, decimal_str, 10);
    int result = mpz_probab_prime_p(n, REPS);
    printf("%s %d\n", decimal_str, result);
    mpz_clear(n);
}

int main(void) {
    /* Edge cases */
    test_value("0");
    test_value("1");
    test_value("2");
    test_value("3");
    test_value("4");

    /* Small primes */
    test_value("5");
    test_value("7");
    test_value("13");
    test_value("97");
    test_value("1009");
    test_value("1000000007");

    /* Small composites */
    test_value("6");
    test_value("8");
    test_value("9");
    test_value("15");
    test_value("100");
    test_value("1000000006");

    /* Strong pseudoprimes to base 2 */
    test_value("2047");     /* 23 * 89 */
    test_value("3277");     /* 29 * 113 */
    test_value("4033");     /* 37 * 109 */
    test_value("4681");     /* 31 * 151 */
    test_value("8321");     /* 53 * 157 */
    test_value("15841");    /* 7 * 31 * 73 */
    test_value("29341");    /* 13 * 37 * 61 */
    test_value("42799");    /* 127 * 337 */
    test_value("52633");    /* 7 * 73 * 103 */
    test_value("65281");    /* 97 * 673 */

    /* Carmichael numbers */
    test_value("561");      /* 3 * 11 * 17 */
    test_value("1105");     /* 5 * 13 * 17 */
    test_value("1729");     /* 7 * 13 * 19 */
    test_value("2465");     /* 5 * 17 * 29 */
    test_value("2821");     /* 7 * 13 * 31 */
    test_value("6601");     /* 7 * 23 * 41 */
    test_value("8911");     /* 7 * 19 * 67 */

    /* Strong pseudoprimes to bases 2 AND 3 */
    test_value("1373653");  /* 829 * 1657 */
    test_value("1530787");

    /* Large prime: 2^127 - 1 (Mersenne) */
    test_value("170141183460469231731687303715884105727");

    /* Large prime: secp256k1 order */
    test_value("115792089237316195423570985008687907852837564279074904382605163141518161494337");

    /* Negative input */
    test_value("-1000000007");

    printf("DONE\n");
    return 0;
}
