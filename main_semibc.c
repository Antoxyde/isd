#include "utils.h"
#include "libpopcnt.h"
#include "semi_bc.h"
#include "xoshiro256starstar.h"

#include <stdint.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>


int main(void) {

    uint64_t time_sec = 2;
    char *challenge_file = "challenges/LW_1280_1";
    uint64_t l = 20;
    uint64_t c = 32;

    printf("# Running with the following configuration :\n");

#if defined(DEBUG)
    printf("# Mode: debug\n");
#else
    printf("# Mode: release\n");
#endif

#if defined(__AVX512DQ__) && defined(__AVX512F__) && defined(__AVX512VL__)
    printf("# AVX512: enabled\n");
#else
    printf("# AVX512: disabled\n");
#endif

    printf("# Algorithm: Semi-BC\n");
    printf("# l: %lu\n", l);
    printf("# Challenge file: %s\n", challenge_file);
    printf("# Time (s): %lu\n", time_sec);

    //uint64_t seed[4] = {1,2,3,4};
    //xoshiro256starstar_random_set(seed);

    uint32_t n = 1280; // Size of the instance
    mzd_t* G = mzd_init(n/2, n);
    mzd_t* H = mzd_init(n/2, n);

    if (!load_challenge(challenge_file, G, H)) {
        return 1;
    }

    mzd_t* min_cw = semi_bc(G, time_sec, l, c);

    if (!min_cw) {
        printf("failed, leaving.\n");
        return 0;
    }

    mzd_t* Hct = mzd_mul(NULL, H, mzd_transpose(NULL, min_cw), 0);
    printf("# Sanity check: %s\n" , mzd_is_zero(Hct) ? "ok" : "nok");

    mzd_free(min_cw);
    mzd_free(G);
    mzd_free(H);
    mzd_free(Hct);
}
