#include "stern.h"
#include "utils.h"
#include "xoshiro256starstar.h"
#include "libpopcnt.h"
#include <time.h>
#include <immintrin.h>

#define STERN_GET(tab, index) ((tab[index >> 6]) >> (index & 0x3f)) & 1
#define STERN_SET_ONE(tab, index) (tab[index >> 6]) |= (1ULL << (index & 0x3f))

#define N 1280
#define K 640

mzd_t* isd_stern_canteaut_chabaud_p2_sort(mzd_t* G, uint64_t time_sec, uint64_t sigma, uint64_t radix_width, uint64_t radix_nlen, uint64_t m, uint64_t c, uint64_t discard_threshold, uint64_t discard_nwords) {

    // Time mesuring stuff
    clock_t start = clock(), current;
    double elapsed = 0.0;

    // p is the Stern p parameter. (= number of LC to make for each rows subsets)
    uint64_t p = 2, iter = 0, mwin = 0, delta = 0, citer = 0, total = 0, done = 0;

    rci_t comb1[p], comb2[p], min_comb[2*p],lambda = 0, mu = 0, tmp = 0, i = 0, j = 0;
    int min_wt = K - 1, wt = 0;

#if defined(AVX512_ENABLED)
    uint64_t mask[10] /* K/64 */;
    memset(mask, 0, 80 /* K/8 */);
#endif

#if defined(DEBUG)
    uint64_t nb_collisions = 0, nb_collisions_delta = 0, collisions_list_size = 0;
#endif

    uint64_t* word = NULL;
    uint64_t* linear_comb = (uint64_t*)malloc(sizeof(uint64_t) * 10 /* K/64 */);
    uint64_t* linear_comb_next = (uint64_t*)malloc(sizeof(uint64_t) * 10 /* K/64 */);
    uint64_t* min_cw = (uint64_t*)malloc(sizeof(uint64_t) * 10 /* K/64 */);
    rci_t* column_perms_copy =  (rci_t*) malloc(sizeof(rci_t) * N);
    rci_t* column_perms = (rci_t*) malloc(sizeof(rci_t) * N);
    CHECK_MALLOC(linear_comb);
    CHECK_MALLOC(linear_comb_next);
    CHECK_MALLOC(min_cw);
    CHECK_MALLOC(column_perms);
    CHECK_MALLOC(column_perms_copy);

    for (i = 0; i < N; i++) column_perms[i] = i;

    mzd_t* Gtemp = mzd_copy(NULL, G);

    // Ensure that we work with a systematic generator matrix
    rref_to_systematic(Gtemp, column_perms);
    mzd_t* Glw = mzd_submatrix(NULL, Gtemp, 0, K, K, N);

    // nelem is the max number of LCs we will make
    // so its an upper bound for the size of each array
    uint64_t nelem = ((N/4 * (N/4 - 1)) /2);

    // Number of bits to store 1 for every possible different values for delta
    uint64_t nb_keys_bits = 1ULL << (sigma - 3);

    // Used to track which combination are present
    uint64_t** collisions_first_pass = (uint64_t**)malloc( m * sizeof(uint64_t*));
    uint64_t** collisions_second_pass = (uint64_t**)malloc( m * sizeof(uint64_t*));

    // Holds the linear combinations made for each pass
    lc** lc_tab_second = (lc**)malloc(m * sizeof(lc*));
    lc** lc_tab_third = (lc**)malloc(m * sizeof(lc*));
    lc** lc_tab_second_sorted = (lc**)malloc(m * sizeof(lc*));
    lc** lc_tab_third_sorted = (lc**)malloc(m * sizeof(lc*));

    CHECK_MALLOC(collisions_first_pass);
    CHECK_MALLOC(collisions_second_pass);
    CHECK_MALLOC(lc_tab_second);
    CHECK_MALLOC(lc_tab_third);
    CHECK_MALLOC(lc_tab_second_sorted);
    CHECK_MALLOC(lc_tab_third_sorted);

    for (mwin = 0; mwin < m; mwin++) {
        collisions_first_pass[mwin] = (uint64_t*)malloc(nb_keys_bits);
        collisions_second_pass[mwin] = (uint64_t*)malloc(nb_keys_bits);
        lc_tab_second[mwin] = (lc*)malloc(nelem * sizeof(lc));
        lc_tab_third[mwin] = (lc*)malloc(nelem * sizeof(lc));
        lc_tab_second_sorted[mwin] = (lc*)malloc(nelem * sizeof(lc));
        lc_tab_third_sorted[mwin] = (lc*)malloc(nelem * sizeof(lc));

        CHECK_MALLOC(collisions_first_pass[mwin]);
        CHECK_MALLOC(collisions_second_pass[mwin]);
        CHECK_MALLOC(lc_tab_second[mwin]);
        CHECK_MALLOC(lc_tab_third[mwin]);
        CHECK_MALLOC(lc_tab_second_sorted[mwin]);
        CHECK_MALLOC(lc_tab_third_sorted[mwin]);
    }

    // These arrays hold the size of their corresponding lc_tab
    uint64_t* lc_tab_second_size = (uint64_t*) malloc(m * sizeof(uint64_t));
    uint64_t* lc_tab_third_size = (uint64_t*) malloc(m * sizeof(uint64_t));
    CHECK_MALLOC(lc_tab_second_size);
    CHECK_MALLOC(lc_tab_third_size);

    // A generic holder for the current indexes
    uint64_t* lc_indexes = (uint64_t*)malloc(m * sizeof(uint64_t));
    CHECK_MALLOC(lc_indexes);

    // Generic alias on the sorted arrays
    lc** lc_tab_alias_second_sorted = (lc**)malloc(m * sizeof(lc*));
    lc** lc_tab_alias_third_sorted = (lc**)malloc(m * sizeof(lc*));

    // Precomputed mask for the window on which we want collision
    uint64_t sigmask = (1ULL << sigma) - 1;

    // Radixsort offsets array
    uint32_t *aux = (uint32_t*) malloc((1ULL << radix_width) * sizeof(uint32_t));

    while (1) {

        for (citer = 0; citer < c; citer++) {

            /* Start of the Canteaut-Chabaud stuff, to derive a new information set. */

            // Find lambda, mu s.t. Glw[lambda, mu] == 1
            uint64_t randdata =  xoshiro256starstar_random();
            lambda = randdata % K;
            randdata >>= 10 /* ceil(log_2(K)) for K=640 */;

            word = Glw->rows[lambda];

            mu = randdata % 10 /* K/64 */;
            randdata >>= 4 /* ceil(log_2(10)) */;

            // Assuming a whole row can't be zero
            while (word[mu] == 0) {
                mu = (mu + 1) % 10 /* K/64 */;
            }

            // We need to take a random 1 in this word.
            // Rotate the word by a random value and take the lowest one it has.
            j  = randdata % 64;
            uint64_t val = ((word[mu] << (64 - j)) | (word[mu] >> j));
            j = (j + _tzcnt_u64(val)) % 64;

#if defined(AVX512_ENABLED)
            uint64_t big_mu = mu, small_mu = j;
#endif
            mu = mu * 64 + j;

           // Log the column swapping
            tmp = column_perms[lambda];
            column_perms[lambda] = column_perms[mu + K];
            column_perms[mu + K] = tmp;

            // Clear the bit at the intersection of the lambda'th row and the mu'th column
            // so we don't have to rewrite a 1 in mu'th column everytime we add the lambda'th row
            mzd_write_bit(Glw, lambda, mu, 0);

#if defined(AVX512_ENABLED)
            void* row_lambda = Glw->rows[lambda];

            // A row contains K bits, so we use 512 + 128 bits registers
            __m512i rlambda1 = _mm512_loadu_si512(row_lambda);
            __m128i rlambda2 = _mm_loadu_si128(row_lambda + 64 /* 512/(8 * sizeof(void)) */);

            mask[big_mu] = (1ULL << small_mu);
            __m512i m1 = _mm512_loadu_si512(mask);
            __m128i m2 = _mm_loadu_si128(((void*)mask) + 64 /* 512/(8 * sizeof(void)) */);
#endif

            // Add the lambda'th row to every other row that have a 1 in the mu'th column
            for (j = 0; j < K; j++) {


#if defined(AVX512_ENABLED)
                void* row = Glw->rows[j];

                // Load the whole row
                __m512i rj1 = _mm512_loadu_si512(row);
                __m128i rj2 = _mm_loadu_si128(row + 64 /* 512/(8 * sizeof(void))  */);

                // Check whether there is a one in column mu using the mask
                if (j != lambda && (_mm512_test_epi64_mask(rj1, m1) >= 1 || _mm_test_epi64_mask(rj2, m2) >= 1)) {

                    // Perform row addition
                    _mm512_storeu_si512(row, _mm512_xor_si512(rlambda1, rj1));
                    _mm_storeu_si128(row + 64 /* 512/(8 * sizeof(void)) */, _mm_xor_si128(rlambda2, rj2));
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
            mask[big_mu] = 0;
#endif
        }

        /* End of the Canteaut-Chabaud stuff. We now have a proper Iset to work with. */

        // Reset all the stuff we will need in the iteration
        memset(lc_indexes, 0, m * sizeof(uint64_t));

        for (mwin = 0; mwin < m; mwin++) {
            memset(collisions_first_pass[mwin], 0, nb_keys_bits);
            memset(collisions_second_pass[mwin], 0, nb_keys_bits);
        }

        // 1st pass, gen all the LC from the first k/2 rows
        for (comb1[0] = 0; comb1[0]  < 320 /* K/2 */; comb1[0]++) {
            uint64_t* row1 = (uint64_t*)Glw->rows[comb1[0]];

            for (comb1[1] = comb1[0] + 1; comb1[1] < 320 /* K/2 */; comb1[1]++) {
                uint64_t* row2 = (uint64_t*)Glw->rows[comb1[1]];

                for (mwin = 0; mwin < m; mwin++) {
                    // Compute the first sigma bits of the LC of rows 1 & 2 on the windows mwin
                    delta = (row1[mwin] ^ row2[mwin]) & sigmask;
                    STERN_SET_ONE(collisions_first_pass[mwin], delta);
                }

            }
        }

        // 2nd pass, gen all the LC from the k/2 to k rows but store
        // only the ones that will collides with at least 1 element from lc_tab_first
        for (comb2[0] = 320 /* K/2 */; comb2[0]  < 640 /* K */; comb2[0]++) {

            uint64_t* row1 = (uint64_t*)Glw->rows[comb2[0]];

            for (comb2[1] = comb2[0] + 1; comb2[1] < 640 /* K */; comb2[1]++) {

                uint64_t* row2 = (uint64_t*)Glw->rows[comb2[1]];

                // The `mwin`st window is the first sigma bits of the `mwin`th word of the row.
                for (mwin = 0; mwin < m; mwin++) {

                    // Compute the first sigma bits of the LC of rows 1 & 2
                    delta = (row1[mwin] ^ row2[mwin]) & sigmask;

                    if (STERN_GET(collisions_first_pass[mwin], delta)) {

                        lc_tab_second[mwin][lc_indexes[mwin]].index1 = comb2[0];
                        lc_tab_second[mwin][lc_indexes[mwin]].index2 = comb2[1];
                        lc_tab_second[mwin][lc_indexes[mwin]].delta = delta;

                        STERN_SET_ONE(collisions_second_pass[mwin], delta);
                        lc_indexes[mwin]++;
                    }
                }
            }
        }

        // Save the size of each lc_tab_first elements and reset lc_indexes
        memcpy(lc_tab_second_size, lc_indexes, m * sizeof(uint64_t));
        memset(lc_indexes, 0, m * sizeof(uint64_t));

        // 3rd pass, copy all the element of first tab that will actually collide in a new array
        for (comb1[0] = 0; comb1[0]  < 320 /* K/2 */; comb1[0]++) {
            uint64_t* row1 = (uint64_t*)Glw->rows[comb1[0]];

            for (comb1[1] = comb1[0] + 1; comb1[1] < 320 /* K/2 */; comb1[1]++) {
                uint64_t* row2 = (uint64_t*)Glw->rows[comb1[1]];

                for (mwin = 0; mwin < m; mwin++) {

                    delta = (row1[mwin] ^ row2[mwin]) & sigmask;

                    if (STERN_GET(collisions_second_pass[mwin], delta)) {
                        lc_tab_third[mwin][lc_indexes[mwin]].index1 = comb1[0];
                        lc_tab_third[mwin][lc_indexes[mwin]].index2 = comb1[1];
                        lc_tab_third[mwin][lc_indexes[mwin]].delta = delta;
                        lc_indexes[mwin]++;
                    }
                }
            }
        }


        // Save the size of each lc_tab_first elements and reset lc_indexes
        memcpy(lc_tab_third_size, lc_indexes, m * sizeof(uint64_t));

        /* From here, lc_tab_third[i][:lc_tab_third_size[i]] and lc_tab_second[i][:lc_tab_second_size[i]]
        *  contains the LCs that WILL collide.
        * So we sort em all. */

        for (mwin = 0; mwin < m; mwin++) {

#if defined(DEBUG)
            collisions_list_size += lc_tab_second_size[mwin] + lc_tab_third_size[mwin];
#endif

            // We have to use aliases to avoid mixing up the lc_tab_X and lc_tab_X_sorted pointers (radixsort returns one of them, depending on the radix_nlen parameter.)
            lc_tab_alias_second_sorted[mwin] = radixsort(lc_tab_second[mwin], lc_tab_second_sorted[mwin], lc_tab_second_size[mwin], radix_width, radix_nlen, aux);
            lc_tab_alias_third_sorted[mwin] = radixsort(lc_tab_third[mwin], lc_tab_third_sorted[mwin],lc_tab_third_size[mwin] , radix_width, radix_nlen, aux);
        }

        // Now that everything is sorted,  we can actually match all the pairs to get the collisions.
        for (mwin = 0; mwin < m; mwin++) {

            int index_second = 0, index_third = 0, save_index_third = 0;

            for (index_second = 0; index_second < lc_tab_second_size[mwin]; index_second++) {

                if (index_second > 0 && lc_tab_alias_second_sorted[mwin][index_second - 1].delta != lc_tab_alias_second_sorted[mwin][index_second].delta) {
                    save_index_third = index_third;

#if defined(DEBUG)
                    nb_collisions_delta++;
#endif
                }

                // Xor the first 2 rows
#if defined(AVX512_ENABLED)
                void* row1 = (void*)Glw->rows[lc_tab_alias_second_sorted[mwin][index_second].index1];
                void* row2 = (void*)Glw->rows[lc_tab_alias_second_sorted[mwin][index_second].index2];

                __m512i linear_comb_high = _mm512_loadu_si512(row1);
                __m128i linear_comb_low = _mm_loadu_si128(row1 + 64 /* 512/8 */ );

                linear_comb_high = _mm512_xor_si512(linear_comb_high, _mm512_loadu_si512(row2));
                linear_comb_low = _mm_xor_si128(linear_comb_low, _mm_loadu_si128(row2 + 64 /* 512/8 */));
#else
                memset(linear_comb, 0, 80 /* K/8 */);
                mxor(linear_comb, (uint64_t*)Glw->rows[lc_tab_alias_second_sorted[mwin][index_second].index1], 10 /* K/64 */);
                mxor(linear_comb, (uint64_t*)Glw->rows[lc_tab_alias_second_sorted[mwin][index_second].index2], 10 /* K/64 */);
#endif

                for (index_third = save_index_third; index_third < lc_tab_third_size[mwin] && lc_tab_alias_second_sorted[mwin][index_second].delta == lc_tab_alias_third_sorted[mwin][index_third].delta; index_third++) {

#if defined(DEBUG)
                    nb_collisions++;
#endif

                    // Xor the two other rows and xor it to the first result
#if defined(AVX512_ENABLED)
                    void* row3 = (void*)Glw->rows[lc_tab_alias_third_sorted[mwin][index_third].index1];
                    void* row4 = (void*)Glw->rows[lc_tab_alias_third_sorted[mwin][index_third].index2];

                    // Load the two new rows and add them to the LC of the two previous ones.
                    __m512i linear_comb_high_next = _mm512_xor_si512(linear_comb_high, _mm512_loadu_si512(row3));
                    __m128i linear_comb_low_next = _mm_xor_si128(linear_comb_low, _mm_loadu_si128(row3 + 64 /* 512/(8 * sizeof(void)) */));

                    linear_comb_high_next = _mm512_xor_si512(linear_comb_high_next, _mm512_loadu_si512(row4));
                    linear_comb_low_next = _mm_xor_si128(linear_comb_low_next, _mm_loadu_si128(row4 + 64 /* 512/(8 * sizeof(void)) */));

                    // Save the result of the LC of the 4 rows
                    _mm512_storeu_si512(linear_comb_next, linear_comb_high_next);
                    _mm_storeu_si128((__m128i*)(linear_comb_next + 8 /* 512/(8 * sizeof(uint64_t) */), linear_comb_low_next);

#else
                    memcpy(linear_comb_next, linear_comb, 80 /* K/8 */ );
                    mxor(linear_comb_next, (uint64_t*)Glw->rows[lc_tab_alias_third_sorted[mwin][index_third].index1], 10 /* K/64 */);
                    mxor(linear_comb_next, (uint64_t*)Glw->rows[lc_tab_alias_third_sorted[mwin][index_third].index2], 10 /* K/64 */);
#endif

                    //printf("DBG mwin = %lu, linear comb is : \n", mwin);
                    //printbin(linear_comb_next, 640);

                    wt = popcnt64_unrolled(linear_comb_next, discard_nwords);
                    total++;

                    // Early abort
                    // For nwords = 7, ~233 is discarding 90% of the codewords
                    if (wt < discard_threshold) {
                        done++;

                        wt += popcnt64_unrolled(linear_comb_next + discard_nwords , 10 /* K/64 */ - discard_nwords);


                        if (wt < min_wt) {

                            // Save the new min weight and the indexes of th e linear combination to obtain it
                            current = clock();
                            elapsed = ((double)(current - start))/CLOCKS_PER_SEC;
                            printf("niter=%lu, time=%.3f, wt=%ld\n", iter, elapsed, wt + 2*p);

                            min_wt = wt;

                            // Save the indexes of the LC
                            min_comb[0] = lc_tab_alias_second_sorted[mwin][index_second].index1;
                            min_comb[1] = lc_tab_alias_second_sorted[mwin][index_second].index2;
                            min_comb[2] = lc_tab_alias_third_sorted[mwin][index_third].index1;
                            min_comb[3] = lc_tab_alias_third_sorted[mwin][index_third].index2;

                            memcpy(min_cw, linear_comb_next, 80 /* K/8 */);
                            memcpy(column_perms_copy, column_perms, N * sizeof(rci_t));

                            mzd_t* cw = stern_reconstruct_cw(min_comb, column_perms_copy, min_cw, p);
                            print_cw(cw);
                            mzd_free(cw);

                            fflush(stdout);
                        }
                    }
                }
            }
        }

        iter++;
        current = clock();
        elapsed = ((double)(current - start))/CLOCKS_PER_SEC;

        if (elapsed > time_sec) {
            break;
        }

    }

#ifdef DEBUG
    printf("# Average number of collisions / iter : %.3f\n", (double)nb_collisions/((double)iter * m));
    printf("# Average number of delta with at least 1 collision / nb delta : %.3f / %lu\n", (double)nb_collisions_delta/((double)iter * m), nelem);
    printf("# Average sorted list size : %.3f\n", (double)collisions_list_size/(2.0 * (double)iter * m));
#endif
    printf("# Total number of iterations done : %lu\n", iter);
    printf("# Iter/s : %.3f\n", (double)iter/(double)time_sec);
    printf("# word analysed, total words, ratio : %lu, %lu, %.3f\n", done, total, (double)done/(double)total);

    mzd_t* result =  stern_reconstruct_cw(min_comb, column_perms_copy, min_cw, p);

    mzd_free(Gtemp);
    mzd_free(Glw);

    free(min_cw);
    free(linear_comb);
    free(linear_comb_next);
    free(column_perms_copy);
    free(column_perms);

    for (mwin = 0; mwin < m; mwin++) {
        free(collisions_first_pass[mwin]);
        free(collisions_second_pass[mwin]);
        free(lc_tab_second[mwin]);
        free(lc_tab_third[mwin]);
        free(lc_tab_second_sorted[mwin]);
        free(lc_tab_third_sorted[mwin]);
    }

    free(lc_tab_second);
    free(lc_tab_third);
    free(lc_tab_second_sorted);
    free(lc_tab_third_sorted);

    free(collisions_first_pass);
    free(collisions_second_pass);

    free(lc_tab_second_size);
    free(lc_tab_third_size);
    free(lc_tab_alias_second_sorted);
    free(lc_tab_alias_third_sorted);

    free(lc_indexes);

    free(aux);

    return result;
}


lc* denomsort_r(lc* T, lc* Ts, int64_t Tlen, uint64_t width, uint64_t pos, uint32_t* Aux) {

    uint32_t mask, k;
    int64_t i;
    k = 1UL << width;
    mask = k - 1;

    memset(Aux, 0, k * sizeof(uint32_t));

    for (i = 0; i < Tlen; i++) {
        Aux[ (T[i].delta >> pos) & mask]++;
    }

    for (i = 1; i < k; i++) {
        Aux[i] += Aux[i - 1];
    }

    for (i = Tlen - 1; i >= 0; i--) {
        uint32_t val = (T[i].delta >> pos) & mask;
        Aux[val]--; Ts[ Aux[val]] = T[i];
    }

    return Ts;
}

lc* radixsort(lc* T, lc* Ts, int64_t Tlen, uint64_t width, uint64_t nlen, uint32_t* aux) {

    int i;
    lc* tmp;

    for (i = 0; i < nlen; i++) {
        tmp = T;
        T = denomsort_r(T, Ts, Tlen, width, i*width, aux);
        Ts = tmp;
    }

    return T;
}
