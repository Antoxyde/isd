
#include "isd.h"
#include "utils.h"
#include "libpopcnt.h"

#include  "xoshiro256starstar.h"


mzd_t* isd_prange_canteaut(mzd_t* G, int niter) {

    rci_t n = G->ncols, i, j,
          row_min_cw = 0,  // retains the row from which the lowest codeword has been found
          lambda, mu, tmp;

    void* row = NULL;
    uint64_t* word = NULL;
    long int wt = 0;

    rci_t* column_perms_copy =  (rci_t*) malloc(sizeof(rci_t) * n);
    rci_t* column_perms = (rci_t*) malloc(sizeof(rci_t) * n);

    if (!column_perms || !column_perms_copy) {
        fprintf(stderr, "Error in %s:  malloc failed.\n", __func__);
        return NULL;
    }
    for (i = 0; i < n; i++) column_perms[i] = i;

    // init to the weight of the 1 vector
    int min_wt = n;
    mzd_t* min_cw = mzd_init(1,n/2);
    mzd_t* Gtemp = mzd_copy(NULL, G);

    // Ensure that we work with a systematic generator matrix
    rref_to_systematic(Gtemp, column_perms);

    // Glw contains only the redundant part of G
    mzd_t* Glw = mzd_submatrix(NULL, Gtemp, 0, n/2, n/2, n);

    for (i = 0; i < niter; i++) {

        // Find lambda, mu s.t. G[lambda, mu] == 1
        lambda = xoshiro256starstar_random() % (n/2);
        do {
            mu = xoshiro256starstar_random() % 10;
            word = mzd_row(Glw, lambda);
        } while (word[mu] == 0);

        word += mu;

        do {
            j  = xoshiro256starstar_random() % 64;
            for (; j < 64 && ((*word >> j) & 1) == 0; j++);
        } while (((*word >> j) & 1) == 0);

        mu = mu*64 + j;

        // Log the column swapping
        tmp = column_perms[lambda];
        column_perms[lambda] = column_perms[mu + (n/2)];
        column_perms[mu + (n/2)] = tmp;

        // Clear the bit at the intersection of the row lambda and the column mu
        // so we don't have to rewrite a 1 in col mu everytime we add the lambda'th row
        mzd_write_bit(Glw, lambda, mu, 0);

        for (j = 0; j < n/2; j++) {
            if (j != lambda && mzd_read_bit(Glw, j, mu) == 1) {

                mzd_row_add(Glw, lambda, j);

                row = mzd_row(Glw, j);
                wt = 1 + popcnt64_unrolled(row, 10 /* 640/64 */ );

                // wt = popcnt(row + (n/16) /* n/2 bits, so n/16 bytes */ , n/16 + (n % 8 != 0) );

                if (wt < min_wt) {
                    printf("New min wt : %ld\n", wt);
                    min_wt = wt;

                    // Save our new lowest row and all the permutations made until now
                    row_min_cw = j;
                    mzd_copy_row(min_cw, 0, Glw, j);
                    memcpy(column_perms_copy, column_perms, n * sizeof(rci_t));
                }
            }
        }

        // Unclear the bit we removed earlier
        mzd_write_bit(Glw, lambda, mu, 1);

    }

    // Since Glw contains only the redundant part of the matrix for the computation
    // we have to concat the identity part back again to get the correct codeword
    mzd_t* ident = mzd_init(1, n/2);

    // row_min_cw contain the row number from which min_cw has been taken
    mzd_write_bit(ident, 0, row_min_cw, 1);
    mzd_t* min_cw_full = mzd_concat(NULL, ident, min_cw);

    // Since we applied many columns permutations, which were all logged in perms
    // we have to permute back to get a valid codeword
    mzd_t* result = mzd_copy(NULL, min_cw_full);
    for (i = 0; i < n; i++) {
        if (i != column_perms_copy[i]) {
            mzd_write_bit(result, 0, column_perms_copy[i], mzd_read_bit(min_cw_full, 0, i));
        }
    }

    mzd_free(Gtemp);
    mzd_free(Glw);
    mzd_free(min_cw);
    mzd_free(min_cw_full);
    mzd_free(ident);

    free(column_perms);
    free(column_perms_copy);

    return result;
}
