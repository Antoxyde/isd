#include "utils.h"
#include "libpopcnt.h"
#include "stern.h"
#include "xoshiro256starstar.h"

#include <stdint.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>


int main(int argc, char** argv) {

    uint64_t niter = 900 * 3600 * 2;
    char *challenge_file = "challenges/LW_1280_1";
    uint64_t sigma = 18;
    uint64_t radix_width = 9;
    uint64_t radix_nlen = 2;

    if (argc > 1) {

        if (argc == 6) {
            challenge_file = argv[1];
            niter = atoi(argv[2]);
            sigma = atoi(argv[3]);
            radix_width = atoi(argv[4]);
            radix_nlen = atoi(argv[5]);
        } else {
            printf("Usage : %s <challenge file> <niter> <sigma> <radix width> <radix nlen>\n", argv[0]);
            return 0;
        }
    }

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

    printf("# Algorithm: Stern\n");
    printf("# Sigma: %lu\n", sigma);
    printf("# Radix width: %lu\n", radix_width);
    printf("# Radix nlen: %lu\n", radix_nlen);
    printf("# Challenge file: %s\n", challenge_file);
    printf("# Niter: %lu\n", niter);

    //uint64_t seed[4] = {1,2,3,4};
    //xoshiro256starstar_random_set(seed);


    clock_t start, stop;
    double time_elapsed;

    uint32_t n = 1280; // Size of the instance
    mzd_t* G = mzd_init(n/2, n);
    mzd_t* H = mzd_init(n/2, n);

    if (!load_challenge(challenge_file, G, H)) {
        return 1;
    }

    start = clock();

    mzd_t* min_cw = isd_stern_canteaut_chabaud_p2_sort(G, niter, sigma, radix_width, radix_nlen);

    if (!min_cw) {
        printf("failed, leaving.\n");
        return 0;
    }

    stop = clock();
    time_elapsed = ((double)(stop - start))/CLOCKS_PER_SEC;

    mzd_t* Hct = mzd_mul(NULL, H, mzd_transpose(NULL, min_cw), 0);

    printf("# Weight: %ld\n", popcnt(mzd_first_row(min_cw), n/8 + (n % 8 != 0)));
    printf("# Verification: %s\n" , mzd_is_zero(Hct) ? "ok" : "nok");
    printf("# Total running time: %.3fs\n", time_elapsed);
    printf("# Iter/s: %.3f\n", ((double)niter)/time_elapsed);

    print_cw(min_cw);
    mzd_free(min_cw);
    mzd_free(G);
    mzd_free(H);
    mzd_free(Hct);
}
