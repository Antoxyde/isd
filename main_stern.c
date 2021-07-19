#include "utils.h"
#include "libpopcnt.h"
#include "stern.h"
#include "xoshiro256starstar.h"

#include <stdint.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

int main(void) {

    uint64_t time_sec = 10;
    char *challenge_file = "challenges/LW_1280_1";

    printf("# Running with the following configuration :\n");

#if defined(DEBUG)
    printf("# Mode              debug\n");
#else
    printf("# Mode              release\n");
#endif

#if defined(__AVX512DQ__) && defined(__AVX512F__) && defined(__AVX512VL__)
    printf("# AVX512            enabled\n");
#else
    printf("# AVX512            disabled\n");
#endif

    printf("# Algorithm         Stern\n");
    printf("# L                 %d\n", L);
    printf("# M                 %d\n", M);
    printf("# P1                %d\n", P1);
    printf("# P2                %d\n", P2);
    printf("# Radixsort width   %d\n", RADIX_WIDTH);
    printf("# Radixsort nlen    %d\n", RADIX_LEN);
    printf("# Challenge file    %s\n", challenge_file);
    printf("# Time (s)          %lu\n", time_sec);

    mzd_t* G = mzd_init(K, N);
    mzd_t* H = mzd_init(K, N);

    if (!load_challenge(challenge_file, G, H)) {
        return 1;
    }

    mzd_t* min_cw = stern(G, time_sec);

    if (!min_cw) {
        printf("failed, leaving.\n");
        return 0;
    }

    mzd_t* Hct = mzd_mul(NULL, H, mzd_transpose(NULL, min_cw), 0);
    printf("# Sanity check      %s\n" , mzd_is_zero(Hct) ? "ok" : "nok");

    mzd_free(min_cw);
    mzd_free(G);
    mzd_free(H);
    mzd_free(Hct);
}
