/*
 * Cross-validation harness for sqisign-verify types, serialization, and hash.
 *
 * Validates: public_key_to_bytes, public_key_from_bytes, signature_to_bytes,
 *   signature_from_bytes, hash_to_challenge.
 *
 * Build: tools/c-validate/build_verify_types.sh
 * Run:   tools/c-validate/verify_types_cval
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

/* Pull in field arithmetic */
#include "fp_p5248_64.c"
#include "fp_select.c"
#include "fp2.c"
#include "mp.c"

/* Precomputed constants */
#include "ec_params.c"

/* EC layer */
#include "ec.c"

/* SHAKE256 implementation */
#include "fips202.c"

/* Verification layer */
#include "encode_verification.c"
#include "common.c"

static void print_hex(const char *label, const unsigned char *data, size_t len)
{
    printf("%s = ", label);
    for (size_t i = 0; i < len; i++)
        printf("%02x", data[i]);
    printf("\n");
}

static void print_fp2_hex(const char *label, const fp2_t *a)
{
    uint8_t buf[FP2_ENCODED_BYTES];
    fp2_encode(buf, a);
    printf("%s = ", label);
    for (int i = 0; i < FP2_ENCODED_BYTES; i++)
        printf("%02x", buf[i]);
    printf("\n");
}

static void print_scalar_hex(const char *label, const scalar_t s)
{
    printf("%s = ", label);
    for (int i = 0; i < NWORDS_ORDER; i++) {
        uint8_t bytes[8];
        memcpy(bytes, &s[i], 8);
        for (int j = 0; j < 8; j++)
            printf("%02x", bytes[j]);
    }
    printf("\n");
}

int main(void)
{
    printf("=== Verify Types Cross-Validation ===\n\n");

    /* --- Section 1: Public key serialization round-trip --- */
    printf("--- Section 1: public_key serialization ---\n");
    {
        public_key_t pk;
        public_key_init(&pk);

        /* Set A = 6 (a simple value), C = 1 (from init) */
        fp2_set_small(&pk.curve.A, 6);
        pk.hint_pk = 0x42;

        unsigned char enc[PUBLICKEY_BYTES];
        public_key_to_bytes(enc, &pk);
        print_hex("pk_bytes", enc, PUBLICKEY_BYTES);

        /* Round-trip */
        public_key_t pk2;
        public_key_from_bytes(&pk2, enc);
        unsigned char enc2[PUBLICKEY_BYTES];
        public_key_to_bytes(enc2, &pk2);
        print_hex("pk_rt_bytes", enc2, PUBLICKEY_BYTES);

        printf("pk_roundtrip = %d\n", memcmp(enc, enc2, PUBLICKEY_BYTES) == 0);
    }

    /* --- Section 2: Signature serialization --- */
    printf("\n--- Section 2: signature_serialization ---\n");
    {
        signature_t sig;
        memset(&sig, 0, sizeof(sig));

        /* Fill with small known values */
        fp2_set_small(&sig.E_aux_A, 42);
        sig.backtracking = 3;
        sig.two_resp_length = 7;

        /* Matrix entries: simple digit values */
        sig.mat_Bchall_can_to_B_chall[0][0][0] = 0x0102030405060708ULL;
        sig.mat_Bchall_can_to_B_chall[0][1][0] = 0x1112131415161718ULL;
        sig.mat_Bchall_can_to_B_chall[1][0][0] = 0x2122232425262728ULL;
        sig.mat_Bchall_can_to_B_chall[1][1][0] = 0x3132333435363738ULL;

        sig.chall_coeff[0] = 0xAABBCCDDEEFF0011ULL;
        sig.hint_aux = 0xAA;
        sig.hint_chall = 0xBB;

        unsigned char enc[SIGNATURE_BYTES];
        signature_to_bytes(enc, &sig);
        print_hex("sig_bytes", enc, SIGNATURE_BYTES);

        /* Round-trip */
        signature_t sig2;
        signature_from_bytes(&sig2, enc);
        unsigned char enc2[SIGNATURE_BYTES];
        signature_to_bytes(enc2, &sig2);
        print_hex("sig_rt_bytes", enc2, SIGNATURE_BYTES);

        printf("sig_roundtrip = %d\n", memcmp(enc, enc2, SIGNATURE_BYTES) == 0);
    }

    /* --- Section 3: hash_to_challenge --- */
    printf("\n--- Section 3: hash_to_challenge ---\n");
    {
        public_key_t pk;
        public_key_init(&pk);
        fp2_set_small(&pk.curve.A, 6);
        pk.hint_pk = 0;

        ec_curve_t com;
        ec_curve_init(&com);
        fp2_set_small(&com.A, 10);

        const unsigned char msg[] = "SQIsign cross-validation test";

        /* Print j-invariants for debugging */
        fp2_t j1, j2;
        ec_j_inv(&j1, &pk.curve);
        ec_j_inv(&j2, &com);
        print_fp2_hex("j1", &j1);
        print_fp2_hex("j2", &j2);

        scalar_t challenge;
        hash_to_challenge(&challenge, &pk, &com, msg, sizeof(msg) - 1);
        print_scalar_hex("challenge", challenge);
    }

    /* --- Section 4: hash_to_challenge with different inputs --- */
    printf("\n--- Section 4: hash_to_challenge variant ---\n");
    {
        public_key_t pk;
        public_key_init(&pk);
        fp2_set_small(&pk.curve.A, 0);
        pk.hint_pk = 0;

        ec_curve_t com;
        ec_curve_init(&com);
        fp2_set_small(&com.A, 0);

        const unsigned char msg[] = "";

        scalar_t challenge;
        hash_to_challenge(&challenge, &pk, &com, msg, 0);
        print_scalar_hex("challenge_zero", challenge);
    }

    printf("\n=== Done ===\n");
    return 0;
}
