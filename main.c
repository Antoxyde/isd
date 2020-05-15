#include "utils.h"
#include "simd.h"
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

int main(void) {

    srand(time(NULL));

    uint32_t n = 1280; // Size of the instance
    mzd_t* G = mzd_init(n/2, n);
    mzd_t* H = mzd_init(n/2, n);

    if (!load_challenge("challenges/LW_1280_0", G, H)) {
        return 1;
    }

    int niter = 100;

    mzd_t* min_cw = isd_lee_brickell(G, niter);

    printf("Min codeword found : \n");

#if defined(__AVX512F__) && defined(__AVX512BW__)
    printf("wt : %d\n", popcnt1280(mzd_first_row(min_cw)));
#else
    printf("wt : %d\n", popcnt(mzd_first_row(min_cw), n/8 + (n % 8 != 0)));
#endif

    mzd_t* Hct = mzd_mul(NULL, H, mzd_transpose(NULL, min_cw), 0);
    printf("Verif : %s\n" , mzd_is_zero(Hct) ? "ok" : "nok");

    mzd_free(min_cw);
    mzd_free(G);
    mzd_free(H);

}
