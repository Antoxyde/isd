#include "isd.h"
#include "utils.h"
#include "libpopcnt.h"

#include <stdint.h>

void get_random_iset(const mzd_t* Gt, mzd_t* Gis, mzd_t* Gist, rci_t* indices) {

    rci_t n = Gt->nrows;

    // M4ri have no function to copy a column from a matrix to another
    // so we transpose and copy rows instead

    // while our matrix is not invertible
    do {

        // we choose n/2 random columns
        fisher_yates_shuffle(indices, n);

        for (rci_t col = 0; col < n/2; col++) {
            mzd_copy_row(Gist, col, Gt, indices[col]);
        }

    } while (mzd_echelonize(Gist, 0) < (n/2 - (n/2) /80) );

    // Since we copied rows from Gt, we have to tranpose back again
    mzd_transpose(Gis, Gist);

}


mzd_t* isd_prange(mzd_t* G, int niter) {

    rci_t n = G->ncols, i = 0;
    mzd_t* Gis = mzd_init(n/2, n/2);
    mzd_t* Gis_inv = mzd_init(n/2, n/2);
    mzd_t* Glw = mzd_init(n/2, n);

    mzd_t* Gt = mzd_transpose(NULL, G);
    mzd_t* Gist = mzd_init(n/2, n/2);

    int min_wt = n;
    mzd_t* min_cw = mzd_init(1,n);

    rci_t* indices = (rci_t*) malloc(sizeof(rci_t) * n);

    if (!indices) {
        fprintf(stderr, "Error in %s: failed to malloc %ld bytes.\n", __func__, sizeof(rci_t) * n);
        return NULL;
    }

    for (i = 0; i < n; i++) indices[i] = i;

    for (i = 0; i < niter; i++) {

        get_random_iset(Gt, Gis, Gist, indices); // Gis is n/2 x n/2

        mzd_inv_m4ri(Gis_inv, Gis, 0);
        mzd_mul(Glw, Gis_inv, G, 0); // Gi * G = Glw  ,  (n/2 x n/2) * (n/2 x n) => n/2 x n

        // Check all the rows of Glw for low codewords
        for (rci_t j = 0; j < n/2; j++) {

            void* row = mzd_row(Glw, j);
            long int wt = popcnt(row, n/8 + (n % 8 != 0) );

            if (wt < min_wt) {
                printf("New min wt : %ld\n", wt);
                min_wt = wt;
                mzd_copy_row(min_cw, 0, Glw, j);
            }
        }

    }

    free(indices);

    mzd_free(Gist);
    mzd_free(Gt);
    mzd_free(Glw);
    mzd_free(Gis);
    mzd_free(Gis_inv);

    return min_cw;
}


mzd_t* isd_lee_brickell(mzd_t* G, int niter) {

    int wt;
    rci_t n = G->ncols, i = 0;
    mzd_t* Gis = mzd_init(n/2, n/2);
    mzd_t* Gis_inv = mzd_init(n/2, n/2);
    mzd_t* Glw = mzd_init(n/2, n);
    mzd_t* Gt = mzd_transpose(NULL, G);
    mzd_t* Gist = mzd_init(n/2, n/2);

    int min_wt = n;
    mzd_t* min_cw = mzd_init(1, n);
    mzd_t* row_k = mzd_init(1, n);
    mzd_t* row_l = mzd_init(1, n);
    mzd_t* row_res = mzd_init(1,n);

    rci_t* indices = (rci_t*) malloc(sizeof(rci_t) * n);

    if (!indices) {
        fprintf(stderr, "Error in %s: failed to malloc %ld bytes.\n", __func__, sizeof(rci_t) * n);
        return NULL;
    }

    for (i = 0; i < n; i++) indices[i] = i;

    for (i = 0; i < niter; i++) {

        get_random_iset(Gt, Gis, Gist, indices); // Gis is n/2 x n/2
        mzd_inv_m4ri(Gis_inv, Gis, 0);
        mzd_mul(Glw, Gis_inv, G, 0); // Gi * G = Glw  ,  (n/2 x n/2) * (n/2 x n) => n/2 x n

        // Check all the CL of the rows for low cw
        for (rci_t k = 0; k < n/2; k++) {

            wt = popcnt(mzd_row(Glw, k), n/8 + (n % 8 != 0) );

            if (wt < min_wt) {
                printf("New min wt : %d\n", wt);
                min_wt = wt;
                mzd_copy_row(min_cw, 0, Glw, k);
            }


            for (rci_t l = 0; l < n/2; l++) {
                if (k != l) {

                    // TODO : optimized CL, without copying to a matrix
                    // eg. using 2 _mm512_loadu_epi8 and one

                    mzd_copy_row(row_l, 0, Glw, l);
                    mzd_copy_row(row_k, 0, Glw, k);
                    mzd_add(row_res, row_k, row_l);

                    wt = popcnt(mzd_first_row(row_res), n/8 + (n % 8 != 0) );

                    if (wt < min_wt) {
                        printf("New min wt : %d\n", wt);
                        min_wt = wt;
                        mzd_copy(min_cw, row_res);
                    }
                }
            }
        }

    }

    free(indices);

    mzd_free(row_l);
    mzd_free(row_k);
    mzd_free(row_res);

    mzd_free(Gt);
    mzd_free(Glw);
    mzd_free(Gist);
    mzd_free(Gis);
    mzd_free(Gis_inv);

    return min_cw;
}


