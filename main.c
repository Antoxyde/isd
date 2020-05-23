#include "utils.h"
#include "libpopcnt.h"
#include "isd.h"

#include <stdint.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

int sanity_check(mzd_t* G, mzd_t* H) {

    int result = 0;
    mzd_t* GHt = mzd_mul(NULL, G, mzd_transpose(NULL, H), 0);
    result += mzd_is_zero(GHt) ? 0 : 1;

    return result == 0;
}

void print_cw(mzd_t* cw) {
    printf("cw = ");
    for (int i = 0; i < cw->ncols; i++) printf(mzd_read_bit(cw, 0, i) ? "1" : "0");
    printf("\n");
}

int main(void) {

    int seed = time(NULL);

    clock_t start, stop;
    double time_elapsed;

    uint32_t n = 1280; // Size of the instance
    int niter = 1;
    mzd_t* G = mzd_init(n/2, n);
    mzd_t* H = mzd_init(n/2, n);

    if (!load_challenge("challenges/LW_1280_0", G, H)) {
        return 1;
    }

    srand(seed);

    start = clock();
    mzd_t* min_cw = isd_prange_canteaut_test(G, niter);
    stop = clock();
    time_elapsed = ((double)(stop - start))/CLOCKS_PER_SEC;

    mzd_t* Hct = mzd_mul(NULL, H, mzd_transpose(NULL, min_cw), 0);

    printf("test version - Min codeword found : \n");
    printf("wt : %ld\n", popcnt(mzd_first_row(min_cw), n/8 + (n % 8 != 0)));
    printf("Verif : %s\n" , mzd_is_zero(Hct) ? "ok" : "nok");
    printf("Total time: %.3f\n", time_elapsed);
    printf("Iter/s : %.3f\n", ((double)niter)/time_elapsed);
    mzd_print(min_cw);

    mzd_free(min_cw);

    srand(seed);

    start = clock();
    min_cw = isd_prange_canteaut(G, niter);
    stop = clock();
    time_elapsed = ((double)(stop - start))/CLOCKS_PER_SEC;

    Hct = mzd_mul(NULL, H, mzd_transpose(NULL, min_cw), 0);

    printf("good version - Min codeword found : \n");
    printf("wt : %ld\n", popcnt(mzd_first_row(min_cw), n/8 + (n % 8 != 0)));
    printf("Verif : %s\n" , mzd_is_zero(Hct) ? "ok" : "nok");
    printf("Total time: %.3f\n", time_elapsed);
    printf("Iter/s : %.3f\n", ((double)niter)/time_elapsed);

    mzd_print(min_cw);

    mzd_free(min_cw);
    mzd_free(G);
    mzd_free(H);
}
