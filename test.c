#include "utils.h"
#include "libpopcnt.h"

#include <stdint.h>
#include <stdlib.h>
#include <time.h>

void fisher_yates_shuffle(rci_t* values, uint32_t size) {

    rci_t i,j, tmp;

    for (i = 0; i < size - 2; i++) {
        j  = rand() % size;
        tmp = values[i];
        values[i] = values[j];
        values[j] = tmp;
    }

}

mzd_t* get_random_iset(mzd_t* G) {

    rci_t n = G->ncols;

    rci_t* indices = (rci_t*) malloc(sizeof(rci_t) * n);

    if (!indices) {
        fprintf(stderr, "Error in %s: failed to malloc %ld bytes.\n", __func__, sizeof(rci_t) * n);
        return NULL;
    }

    for (rci_t i = 0; i < n; i++) indices[i] = i;

    // the goal is to find an invertible n/2 \times n/2 matrix
    mzd_t* G1 = mzd_init(n/2, n/2);

    // M4ri have no function to copy a column from a matrix to another
    // so we transpose and copy rows instead
    mzd_t* Gt = mzd_transpose(NULL, G);

    // while our matrix is not invertible
    do {

        // we choose n/2 random columns
        fisher_yates_shuffle(indices, n);

        for (rci_t col = 0; col < n/2; col++) {
            mzd_copy_row(G1, col, Gt, indices[col]);
        }

    } while (mzd_echelonize(G1, 0) != n/2);

    // Since we copied rows from Gt, we have to tranpose back again
    mzd_t* Gi = mzd_transpose(NULL, G1);

    mzd_free(G1);
    mzd_free(Gt);
    free(indices);

    return Gi;
}


int sanity_check(mzd_t* G, mzd_t* H) {

    int result = 0;
    mzd_t* GHt = mzd_mul(NULL, G, mzd_transpose(NULL, H), 0);
    result += mzd_is_zero(GHt) ? 0 : 1;

    return result == 0;
}

int main(void) {

    srand(time(NULL)); // chosen by fair dice roll

    uint32_t n = 1280; // Size of the instance
    mzd_t* G = mzd_init(n/2, n);
    mzd_t* H = mzd_init(n/2, n);
    load_challenge("challenges/LW_1280_0", G, H);

    mzd_t* Gi = get_random_iset(G); // n/2 x n/2

    mzd_t* Glw = mzd_mul(NULL, Gi, G, 0); // (n/2 x n/2) * (n/2 x n) => n/2 x n

    int min_wt = n;
    mzd_t* min_cw = mzd_init(1,n);

    // Check all the rows of Glw for low codewords
    for (rci_t i = 0; i < Glw->nrows; i++) {
        void* row = mzd_row(Glw, i);
        int wt = popcnt(row, n/8);

        if (wt < min_wt) {
            min_wt = wt;
            mzd_copy_row(min_cw, 0, Glw, i);
        }
    }

    printf("Lowest codeword found :\n");
    mzd_print(min_cw);

    printf("wt : %d\n", min_wt);

    mzd_t* Hct = mzd_mul(NULL, H, mzd_transpose(NULL, min_cw), 0);
    printf("Verif : %s\n" , mzd_is_zero(Hct) ? "ok" : "nok");

    mzd_free(min_cw);
    mzd_free(Hct);
    mzd_free(Glw);
    mzd_free(Gi);
    mzd_free(G);
    mzd_free(H);

    return 0;
}
