#include "utils.h"
#include "libpopcnt.h"
#include <stdint.h>

/*
rci_t* get_random_iset(mzd_t* H) {

}
*/

int main(void) {

    uint32_t n = 1280; // Size of the instance

    mzd_t* G = mzd_init(n/2, n);
    mzd_t* H = mzd_init(n/2, n);

    load_challenge("challenges/LW_1280_0", G, H);

    mzd_t* GHt = mzd_mul_m4rm(NULL, G, mzd_transpose(NULL, H), 0);
    printf("GH^t : \n");
    mzd_print(GHt); // Devrait être NULL

    printf("G : \n");
    mzd_echelonize(G, 1);
    //mzd_print(G);


/*
    uint64_t count = popcnt((void*)mzd_first_row(H), 640/8);
    printf("First half line popcnt: %lu.\n", count); // Should be 1

    uint64_t rank = mzd_echelonize(H, 1);
    printf("Matrix rank: %lu.\n", rank); // Should be 640

*/

    return 0;
}
