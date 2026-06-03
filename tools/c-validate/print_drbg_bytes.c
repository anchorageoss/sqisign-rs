/*
 * Print the CTR-DRBG bytes that biextension-test.c would consume with seed=0.
 */
#include <stdint.h>
#include <stdio.h>
#include <string.h>

#include "aes_c.c"

#define RANDOMBYTES_C
#define SQISIGN_API
#include "randombytes_ctrdrbg.c"

int main(void) {
    uint32_t seed[12] = {0};
    randombytes_init((unsigned char *)seed, NULL, 256);

    /* The test draws 4 x 24 bytes (NWORDS_ORDER-1 = 3 words = 24 bytes) */
    unsigned char buf[24];
    const char *names[] = {"scal_d1", "scal_d2", "scal_s1", "scal_s2"};
    for (int t = 0; t < 4; t++) {
        randombytes(buf, 24);
        printf("%s = ", names[t]);
        for (int i = 0; i < 24; i++) printf("%02x", buf[i]);
        printf("\n");
    }
    return 0;
}
