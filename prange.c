
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
    rci_t k = n/2;

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
    mzd_t* min_cw = mzd_init(1,k);
    mzd_t* Gtemp = mzd_copy(NULL, G);

    // Ensure that we work with a systematic generator matrix
    rref_to_systematic(Gtemp, column_perms);

    // Glw contains only the redundant part of G
    mzd_t* Glw = mzd_submatrix(NULL, Gtemp, 0, k, k, n);

    for (iter = 0; iter < niter; iter++) {

        // Find lambda, mu s.t. Glw[lambda, mu] == 1
        lambda = xoshiro256starstar_random() % 640 /* k */;
        word = Glw->rows[lambda];

        mu = xoshiro256starstar_random() % 10;
        // Assuming a whole row can't be zero
        while (word[mu] == 0) {
            mu = (mu + 1) % 10;
        }

        j  = xoshiro256starstar_random() % 64;
        uint64_t val = ((word[mu] << (64 - j)) | (word[mu] >> j));
        j = (j + _tzcnt_u64(val)) % 64;

#if defined(AVX512_ENABLED)
        uint64_t big_mu = mu, small_mu = j;
#endif
        mu = mu * 64 + j;

        // Log the column swapping
        tmp = column_perms[lambda];
        column_perms[lambda] = column_perms[mu + 640 /* k */];
        column_perms[mu + (n/2)] = tmp;

        // Clear the bit at the intersection of the lambda'th row and the mu'th column
        // so we don't have to rewrite a 1 in mu'th column everytime we add the lambda'th row
        mzd_write_bit(Glw, lambda, mu, 0);

#if defined(AVX512_ENABLED)
        //void* row_lambda = mzd_row(Glw, lambda);
        void* row_lambda = Glw->rows[lambda];
        __m512i rlambda1 = _mm512_loadu_si512(row_lambda);
        __m128i rlambda2 = _mm_loadu_si128(row_lambda + 64 /* = 512 / (8 * sizeof(void)) */);

        mask[big_mu] = (1ULL << small_mu);
        __m512i m1 = _mm512_loadu_si512(mask);
        __m128i m2 = _mm_loadu_si128(((void*)mask) + 64 /* = 512/(8 * sizeof(void)) */);
#endif

        // Add the lambda'th row to every other row that have a 1 in the mu'th column
        for (j = 0; j < 640 /* k */; j++) {
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
                    min_wt = wt;

                    // Save our new lowest row and all the permutations made until now
                    row_min_cw = j;
                    mzd_copy_row(min_cw, 0, Glw, j);
                    memcpy(column_perms_copy, column_perms, 1280 /* n */ * sizeof(rci_t));

                    mzd_t* cw = prange_reconstruct_cw(row_min_cw, column_perms_copy, min_cw);
                    print_cw(cw);
                    mzd_free(cw),
                    fflush(stdout);

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

    mzd_t* result = prange_reconstruct_cw(row_min_cw, column_perms_copy, min_cw);

    mzd_free(Gtemp);
    mzd_free(Glw);
    mzd_free(min_cw);

    free(column_perms);
    free(column_perms_copy);

    return result;
}
