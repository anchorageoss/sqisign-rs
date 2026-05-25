/*
 * Standalone cross-validation harness for sqisign-fp Level 1.
 *
 * Compiles only the static helpers from fp_p5248_64.c (lines 1..518)
 * by directly #including that file with the bottom half excluded via
 * a sentinel macro. Performs a fixed sequence of operations on known
 * inputs and prints each result as 32 little-endian hex bytes. The
 * Rust test crates/fp/tests/c_crossvalidate.rs runs the same
 * sequence in Rust and compares byte for byte.
 *
 * Build:  gcc -O2 -o cval cval.c
 * Run:    ./cval
 */

#include <stdint.h>
#include <stdio.h>
#include <string.h>

#define spint uint64_t
#define sspint int64_t
#define dpint __uint128_t
#define udpint __uint128_t

/* Forward-declare the static helpers we will use below. */
static void modadd(const spint *a, const spint *b, spint *n);
static void modsub(const spint *a, const spint *b, spint *n);
static void modneg(const spint *b, spint *n);
static void modmul(const spint *a, const spint *b, spint *c);
static void modsqr(const spint *a, spint *c);
static void modinv(const spint *x, const spint *h, spint *z);
static void modsqrt(const spint *x, const spint *h, spint *r);
static int modfsb(spint *n);
static void modshl(unsigned int n, spint *a);
static int modshr(unsigned int n, spint *a);
static void nres(const spint *m, spint *n);
static void redc(const spint *n, spint *m);
static void modcpy(const spint *a, spint *c);

/* Include just the first 522 lines of fp_p5248_64.c, the part that
 * defines the static helpers we use, with no external dependencies.
 * The truncated copy is generated at build time by the harness
 * shell script (build.sh). */
#include "fp_p5248_64_static.c"

/* Encode helper: 5-limb internal -> 32 canonical LE bytes. Mirrors
 * fp_encode body (modified modexp).  */
static void encode32(uint8_t *dst, const spint *a)
{
    spint c[5];
    redc(a, c);
    for (int i = 0; i < 32; i++) {
        dst[i] = (uint8_t)(c[0] & 0xff);
        (void)modshr(8, c);
    }
}

/* Decode helper: 32 canonical LE bytes -> 5-limb internal Montgomery
 * form. Mirrors fp_decode body. */
static int decode32(spint *d, const uint8_t *src)
{
    int i;
    for (i = 0; i < 5; i++) d[i] = 0;
    for (i = 31; i >= 0; i--) {
        modshl(8, d);
        d[0] += (spint)src[i];
    }
    int in_range = modfsb(d);
    spint tmp[5];
    nres(d, tmp);
    /* fp_decode masks the output to zero if out of range; for our
     * test we only feed in-range inputs. */
    (void)in_range;
    /* nres writes to `d` from its second arg `tmp`, but the actual
     * C `nres` reads `m` then writes `n`. Re-do correctly: */
    nres(d, d);  /* OK -- nres aliases in/out via modmul */
    return in_range;
}

static void print_hex(const char *label, const uint8_t *bytes)
{
    printf("%s = ", label);
    for (int i = 0; i < 32; i++) {
        printf("%02x", bytes[i]);
    }
    printf("\n");
}

int main(void)
{
    /* Three canonical inputs (LE bytes), chosen to exercise the
     * low/high regions of the prime. */
    uint8_t a_bytes[32] = {
        0x01, 0x23, 0x45, 0x67, 0x89, 0xab, 0xcd, 0xef,
        0xfe, 0xdc, 0xba, 0x98, 0x76, 0x54, 0x32, 0x10,
        0x00, 0x11, 0x22, 0x33, 0x44, 0x55, 0x66, 0x77,
        0x88, 0x99, 0xaa, 0xbb, 0xcc, 0xdd, 0xee, 0x03
    };
    uint8_t b_bytes[32] = {
        0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
        0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
        0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
        0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0x03
    };

    spint a[5], b[5], r[5];
    decode32(a, a_bytes);
    decode32(b, b_bytes);

    uint8_t out[32];

    print_hex("a", a_bytes);
    print_hex("b", b_bytes);

    modadd(a, b, r);
    encode32(out, r);
    print_hex("a+b", out);

    modsub(a, b, r);
    encode32(out, r);
    print_hex("a-b", out);

    modneg(a, r);
    encode32(out, r);
    print_hex("-a", out);

    modmul(a, b, r);
    encode32(out, r);
    print_hex("a*b", out);

    modsqr(a, r);
    encode32(out, r);
    print_hex("a^2", out);

    spint t[5];
    modcpy(a, t);
    modinv(t, NULL, r);
    encode32(out, r);
    print_hex("1/a", out);

    /* For sqrt, compute a^2 first (always a QR), then sqrt back. */
    modsqr(a, t);
    modsqrt(t, NULL, r);
    encode32(out, r);
    print_hex("sqrt(a^2)", out);

    return 0;
}
