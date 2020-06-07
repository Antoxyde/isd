
#include "prange.h"
#include "utils.h"
#include "libpopcnt.h"
#include  "xoshiro256starstar.h"

#include <time.h>


mzd_t* isd_prange_canteaut_chabaud(mzd_t* G, uint64_t niter) {

    clock_t start, current;
    start = clock();

    rci_t n = G->ncols, i, j,
          row_min_cw = 0,  // retains the row from which the lowest codeword has been found
          lambda, mu, tmp;

    void* row = NULL;
    uint64_t* word = NULL;
    uint64_t iter;

#if defined(AVX512_ENABLED)
    uint64_t mask[10];
    memset(mask, 0, 80);
#endif

    long int wt = 0;

    rci_t* column_perms_copy =  (rci_t*) malloc(sizeof(rci_t) * n);
    rci_t* column_perms = (rci_t*) malloc(sizeof(rci_t) * n);
    CHECK_MALLOC(column_perms);
    CHECK_MALLOC(column_perms_copy);
    for (i = 0; i < n; i++) column_perms[i] = i;

    // init to the weight of the 1 vector
    int min_wt = n;
    mzd_t* min_cw = mzd_init(1,n/2);
    mzd_t* Gtemp = mzd_copy(NULL, G);

    // Ensure that we work with a systematic generator matrix
    rref_to_systematic(Gtemp, column_perms);

    // Glw contains only the redundant part of G
    mzd_t* Glw = mzd_submatrix(NULL, Gtemp, 0, n/2, n/2, n);

    for (iter = 0; iter < niter; iter++) {

        // Find lambda, mu s.t. Glw[lambda, mu] == 1
        // Assuming that a whole row can't be totally 0, but that 64 bits subset of that row
        lambda = xoshiro256starstar_random() % (n/2);
        do {
            mu = xoshiro256starstar_random() % 10;
            //word = mzd_row(Glw, lambda);
            word = Glw->rows[lambda];
        } while (word[mu] == 0);

        word += mu;

        do {
            j  = xoshiro256starstar_random() % 64;
            for (; j < 64 && ((*word >> j) & 1) == 0; j++);
        } while (((*word >> j) & 1) == 0);

        mu = mu * 64 + j;

        // Log the column swapping
        tmp = column_perms[lambda];
        column_perms[lambda] = column_perms[mu + (n/2)];
        column_perms[mu + (n/2)] = tmp;

        // Clear the bit at the intersection of the lambda'th row and the mu'th column
        // so we don't have to rewrite a 1 in mu'th column everytime we add the lambda'th row
        mzd_write_bit(Glw, lambda, mu, 0);

#if defined(AVX512_ENABLED)
        //void* row_lambda = mzd_row(Glw, lambda);
        void* row_lambda = Glw->rows[lambda];
        __m512i rlambda1 = _mm512_loadu_si512(row_lambda);
        __m128i rlambda2 = _mm_loadu_si128(row_lambda + 64 /* = 512 / (8 * sizeof(void)) */);

        // No easy instrinsic to set a single bit to 1 ?
        mask[mu/64] = ((uint64_t)1 << (mu%64));
        __m512i m1 = _mm512_loadu_si512(mask);
        __m128i m2 = _mm_loadu_si128(((void*)mask) + 64 /* = 512/(8 * sizeof(void)) */);
#endif

        // Add the lambda'th row to every other row that have a 1 in the mu'th column
        for (j = 0; j < n/2; j++) {
            //row = mzd_row(Glw, j);
            row = Glw->rows[j];

#if defined(AVX512_ENABLED)
            // Load the whole row
            __m512i rj1 = _mm512_loadu_si512(row);
            __m128i rj2 = _mm_loadu_si128(row + 64 /* = 512 / (8 * sizeof(void)) */);

            // Check whether there is a one in column mu using the mask
            if (j != lambda && (_mm512_test_epi64_mask(rj1, m1) >= 1 || _mm_test_epi64_mask(rj2, m2) >= 1)) {

                // Perform row addition
                _mm512_storeu_si512(row, _mm512_xor_si512(rlambda1, rj1));
                _mm_storeu_si128(row + 64, _mm_xor_si128(rlambda2, rj2));
#else
            if (j != lambda && mzd_read_bit(Glw, j, mu) == 1) {
                mzd_row_add(Glw, lambda, j);
#endif

                wt = 1 + popcnt64_unrolled(row, 10 /* 640/64 */ );
                // wt = popcnt(row + (n/16) /* n/2 bits, so n/16 bytes */ , n/16 + (n % 8 != 0) );

                if (wt < min_wt) {
                    current = clock();
                    double elapsed = ((double)(current - start))/CLOCKS_PER_SEC;
                    printf("niter=%lu, time=%.3f, wt=%ld\n", iter, elapsed, wt);
                    fflush(stdout);
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

#if defined(AVX512_ENABLED)
        // Clear the mask for the next iteration
        mask[mu/64] = 0;
#endif

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