#include "isd.h"
#include "utils.h"
#include "libpopcnt.h"

mzd_t* isd_prange(mzd_t* G, int niter) {

    rci_t n = G->ncols, i = 0;
    mzd_t* Gis = mzd_init(n/2, n/2);
    mzd_t* Glw = mzd_init(n/2, n);

    int min_wt = n;
    mzd_t* min_cw = mzd_init(1,n);

    rci_t* indices = (rci_t*) malloc(sizeof(rci_t) * n);

    if (!indices) {
        fprintf(stderr, "Error in %s: failed to malloc %ld bytes.\n", __func__, sizeof(rci_t) * n);
        return NULL;
    }

    for (i = 0; i < n; i++) indices[i] = i;

    for (i = 0; i < niter; i++) {

        get_random_iset(G, Gis, indices); // Gis is n/2 x n/2
        mzd_mul(Glw, Gis, G, 0); // Gi * G = Glw  ,  (n/2 x n/2) * (n/2 x n) => n/2 x n

        // Check all the rows of Glw for low codewords
        for (rci_t j = 0; j < n/2; j++) {
            void* row = mzd_row(Glw, j);
            int wt = popcnt(row, n/8);

            if (wt < min_wt) {
                printf("New min wt : %d\n", wt);
                min_wt = wt;
                mzd_copy_row(min_cw, 0, Glw, j);
            }
        }

    }

    free(indices);

    mzd_free(Glw);
    mzd_free(Gis);

    return min_cw;
}

void get_random_iset(mzd_t* G, mzd_t* Gis, rci_t* indices) {

    rci_t n = G->ncols;

    // M4ri have no function to copy a column from a matrix to another
    // so we transpose and copy rows instead
    mzd_t* Gt = mzd_transpose(NULL, G);
    mzd_t* Gist = mzd_init(n/2, n/2);

    // while our matrix is not invertible
    do {

        // we choose n/2 random columns
        fisher_yates_shuffle(indices, n);

        for (rci_t col = 0; col < n/2; col++) {
            mzd_copy_row(Gist, col, Gt, indices[col]);
        }

    } while (mzd_echelonize(Gist, 0) != n/2);

    // Since we copied rows from Gt, we have to tranpose back again
    mzd_transpose(Gis, Gist);

    mzd_free(Gt);
    mzd_free(Gist);
}




