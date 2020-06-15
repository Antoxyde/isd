#include "utils.h"
#include "libpopcnt.h"
#include "prange.h"
#include "xoshiro256starstar.h"

#include <stdint.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

int main(int argc, char** argv) {


    uint64_t niter = 1000;
    char *challenge_file = "challenges/LW_1280_1";

    if (argc > 1) {

        if (argc == 3) {
            challenge_file = argv[1];
            niter = atoi(argv[2]);
        } else {
            printf("Usage : %s <challenge file> <niter>\n", argv[0]);
            return 0;
        }
    }


    //uint64_t seed[4] = {1,2,3,4};
    //xoshiro256starstar_random_set(seed);

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

    printf("# Algorithm: Prange\n");
    printf("# Challenge file: %s\n", challenge_file);
    printf("# Niter: %lu\n", niter);

    clock_t start, stop;
    double time_elapsed;

    uint32_t n = 1280; // Size of the instance
    mzd_t* G = mzd_init(n/2, n);
    mzd_t* H = mzd_init(n/2, n);

    if (!load_challenge(challenge_file, G, H)) {
        return 1;
    }

    start = clock();
    mzd_t* min_cw = isd_prange_canteaut_chabaud(G, niter);

    if (!min_cw) {
        printf("failed, leaving.\n");
        return 0;
    }

    stop = clock();
    time_elapsed = ((double)(stop - start))/CLOCKS_PER_SEC;

    mzd_t* Hct = mzd_mul(NULL, H, mzd_transpose(NULL, min_cw), 0);

    printf("# Sanity check: %s\n" , mzd_is_zero(Hct) ? "ok" : "nok");
    printf("# Total running time: %.3f\n", time_elapsed);
    printf("# Iter/s : %.3f\n", ((double)niter)/time_elapsed);

    mzd_free(min_cw);
    mzd_free(G);
    mzd_free(H);
    mzd_free(Hct);
}
