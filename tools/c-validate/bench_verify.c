#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include <sqisign_namespace.h>
#include <api.h>

static int hex_to_byte(char c) {
    if (c >= '0' && c <= '9') return c - '0';
    if (c >= 'a' && c <= 'f') return c - 'a' + 10;
    if (c >= 'A' && c <= 'F') return c - 'A' + 10;
    return -1;
}

static int hex_to_bytes(const char *hex, unsigned char *out, size_t max_len) {
    size_t hex_len = strlen(hex);
    size_t byte_len = hex_len / 2;
    if (byte_len > max_len) byte_len = max_len;
    for (size_t i = 0; i < byte_len; i++) {
        int hi = hex_to_byte(hex[2*i]);
        int lo = hex_to_byte(hex[2*i+1]);
        if (hi < 0 || lo < 0) return -1;
        out[i] = (unsigned char)((hi << 4) | lo);
    }
    return (int)byte_len;
}

int main(int argc, char **argv) {
    const char *kat_path = "../../reference/KAT/PQCsignKAT_353_SQIsign_lvl1.rsp";
    if (argc > 1) kat_path = argv[1];

    FILE *f = fopen(kat_path, "r");
    if (!f) {
        fprintf(stderr, "Cannot open %s\n", kat_path);
        return 1;
    }

    unsigned char pk[CRYPTO_PUBLICKEYBYTES];
    unsigned char sm[8192];
    size_t smlen = 0;
    int have_pk = 0, have_sm = 0;

    char line[65536];
    while (fgets(line, sizeof(line), f)) {
        if (strncmp(line, "pk = ", 5) == 0) {
            char *hex = line + 5;
            hex[strcspn(hex, "\r\n")] = 0;
            hex_to_bytes(hex, pk, sizeof(pk));
            have_pk = 1;
        } else if (strncmp(line, "sm = ", 5) == 0) {
            char *hex = line + 5;
            hex[strcspn(hex, "\r\n")] = 0;
            smlen = strlen(hex) / 2;
            hex_to_bytes(hex, sm, sizeof(sm));
            have_sm = 1;
        }
        if (have_pk && have_sm) break;
    }
    fclose(f);

    if (!have_pk || !have_sm) {
        fprintf(stderr, "Failed to parse KAT entry\n");
        return 1;
    }

    unsigned char msg[8192];
    unsigned long long msglen = 0;

    /* Warm-up */
    for (int i = 0; i < 5; i++) {
        int ret = crypto_sign_open(msg, &msglen, sm, smlen, pk);
        if (ret != 0) {
            fprintf(stderr, "Verification failed!\n");
            return 1;
        }
    }

    /* Timed run */
    int iterations = 100;
    struct timespec start, end;
    clock_gettime(CLOCK_MONOTONIC, &start);
    for (int i = 0; i < iterations; i++) {
        crypto_sign_open(msg, &msglen, sm, smlen, pk);
    }
    clock_gettime(CLOCK_MONOTONIC, &end);

    double elapsed = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;
    double per_verify_ms = (elapsed / iterations) * 1000.0;

    printf("C reference (ref, Level 1):\n");
    printf("  %d iterations in %.3f s\n", iterations, elapsed);
    printf("  %.3f ms per verification\n", per_verify_ms);

    return 0;
}
