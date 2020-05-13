#include "utils.h"
#include "libpopcnt.h"
#include "isd.h"

#include <stdint.h>
#include <stdlib.h>
#include <time.h>


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
    load_challenge("challenges/LW_1280_0", G, H);

    mzd_t* min_cw = isd_prange(G, 100);

    printf("Min codeword found : \n");
    mzd_print(min_cw);

    printf("wt : %ld\n", popcnt(mzd_first_row(min_cw), n/8));

    mzd_t* Hct = mzd_mul(NULL, H, mzd_transpose(NULL, min_cw), 0);
    printf("Verif : %s\n" , mzd_is_zero(Hct) ? "ok" : "nok");

    mzd_free(min_cw);
    mzd_free(G);
    mzd_free(H);

}

