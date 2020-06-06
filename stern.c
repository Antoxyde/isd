#include "stern.h"
#include "table.h"
#include "utils.h"

mzd_t* isd_stern_canteaut_chabaud_p2(mzd_t* G, uint64_t niter, uint64_t sigma) {

    uint64_t p = 2;
    rci_t n = G->ncols, comb1[2], comb2[2];

    uint64_t* linear_comb = (uint64_t*)malloc(sizeof(uint64_t) * 10);
    uint64_t* min_cw = (uint64_t*)malloc(sizeof(uint64_t) * 10);
    rci_t* column_perms_copy =  (rci_t*) malloc(sizeof(rci_t) * n);
    rci_t* column_perms = (rci_t*) malloc(sizeof(rci_t) * n);
    CHECK_MALLOC(colmun_perms);
    CHECK_MALLOC(column_perms_copy);
    CHECK_MALLOC(linear_comb);
    CHECK_MALLOC(min_cw);

    for (i = 0; i < n; i++) column_perms[i] = i;

    mzd_t* Gtemp = mzd_copy(NULL, G);

    // Ensure that we work with a systematic generator matrix
    rref_to_systematic(Gtemp, column_perms);
    mzd_t* Glw = mzd_submatrix(NULL, Gtemp, 0, n/2, n/2, n);

    // We store v=(a,b,c),k=(Z[a]+Z[b]+Z[c])[:sigma]
    table* tab = table_init(1 << sigma, 1);

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
            }
        }

         // Unclear the bit we removed earlier
        mzd_write_bit(Glw, lambda, mu, 1);

#if defined(AVX512_ENABLED)
        // Clear the mask for the next iteration
        mask[mu/64] = 0;
#endif

        // Here we have a new iset, we can start a Stern iteration

        uint64_t delta1, delta2;
        int min_wt = n;

        for (comb1[0] = 0; comb1[0]  < n/4; comb1[0]++) {
            for (comb1[1] = 0; comb1[1] < n/4; comb1[1]++) {
                if (comb1[0] != comb1[1]) {
                    delta1 = mxor(mzd_row(Glw, comb1[0]), mzd_row(Glw, comb1[1]), sigma);
                    table_insert(tab, &comb1, p, delta1);
                }
            }
        }

        for (comb2[0] = n/4; comb2[0]  < n/2; comb2[0]++) {
            for (comb2[1] = n/4; comb2[1] < n/2; comb2[1]++) {
                if (comb2[0] != comb2[1]) {

                    // Compute the key of the linear combination
                    delta2 = uxor(mzd_row(Glw, comb2[0]), mzd_row(Glw, comb2[1]), sigma);

                    // And check if some elements from the previous set already had this key
                    bucket* buck = table_retrieve_bucket(ht, delta2);

                    if (buck) {
                        for (i = 0; i < buck->len; i++) {
                            // For each element that had the same key, we have a collision
                            if (buck->elems[i]->key == delta2) {

#if defined(AVX512_ENABLED)
                                // TODO cl av512
#else
                                comb1 = (rci_t*)(buck->elems[i]->data);
                                mxor(linear_comb,linear_comb, 80);
                                mxor(linear_comb,mzd_row(Glw, comb1[0]), 80);
                                mxor(linear_comb, mzd_row(Glw, comb1[1]), 80);
                                mxor(linear_comb, mzd_row(Glw, comb2[0]), 80);
                                mxor(linear_comb, mzd_row(Glw, comb2[1]), 80);
#endif
                                // TODO DBG : check linear_comb[:sigma] est bien a zéro
                                int wt = 2*p + popcnt(linear_comb, 10);

                                if (wt < min_wt) {
                                    // Save the new min weight and the indexes of th e linear combination to obtain it
                                    min_wt = wt;
                                    memcpy(min_comb, comb1, p * sizeof(rci_t));
                                    mempcy(min_comb + p, comb2, p * sizeof(rci_t));
                                    memcpy(min_cw, linear_comb, 80);
                                    memcpy(column_perms_copy, column_perms, n * sizeof(rci_t));
                                }
                            }
                        }
                    }
                }
            }
        }


        // Reset the table for the next iteration
        table_reset(t);

    }

    mzd_t* ident = mzd_init(1, n/2);
    for (int i = 0; i < 2*p; i++)
        mzd_write_bit(ident, 0, min_comb[i], 1);

    mzd_t* min_cw_m = mzd_init(1, n/2);
    memcpy(mzd_first_row(min_cw_m), min_cw, 80);

    mzd_t* min_cw_full = mzd_concat(NULL, ident, min_cw_m);

    mzd_t* result = mzd_copy(NULL, min_cw_full);
    for (i = 0; i < n; i++) {
        if (i != column_perms_copy[i]) {
            mzd_write_bit(result, 0, column_perms_copy[i], mzd_read_bit(min_cw_full, 0, i));
        }
    }

    mzd_free(Gtemp);
    mzd_free(Glw);
    mzd_free(min_cw_full);
    mzd_free(ident);

    free(min_cw);
    free(linear_comb);
    free(column_perms_copy);
    free(column_perms);

    return result;
}
