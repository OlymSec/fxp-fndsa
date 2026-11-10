#define FALCON_MQ_TEST 1
#include "../vendor/falcon/vrfy.c"
#include <assert.h>
#include <stdio.h>
static unsigned state = 1;
static unsigned r32(void)
{
    state ^= state << 13;
    state ^= state >> 17;
    state ^= state << 5;
    return state;
}
int main(void)
{
    unsigned count = 0;
    for (unsigned logn = 1; logn <= 10; logn++)
        for (unsigned k = 0; k < 300; k++) {
            unsigned n = 1U << logn;
            uint16_t a[1024], b[1024], original[1024];
            for (unsigned i = 0; i < n; i++)
                a[i] = b[i] = original[i] = (uint16_t)(r32() % 12289);
            mq_NTT(a, logn);
            mq_NTT_reference(b, logn);
            assert(!memcmp(a, b, n * sizeof *a));
            mq_iNTT(a, logn);
            mq_iNTT_reference(b, logn);
            assert(!memcmp(a, b, n * sizeof *a));
            assert(!memcmp(a, original, n * sizeof *a));
            count++;
        }
    printf("PASS: %u forward/inverse transforms match the scalar arithmetic exactly\n", count);
    return 0;
}
