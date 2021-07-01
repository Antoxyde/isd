#include "semi_bc.h"
#include "utils.h"
#include "xoshiro256starstar.h"
#include "libpopcnt.h"
#include <time.h>
#include <immintrin.h>

#define N 1280
#define K 640
#define RADIX_WIDTH 10
#define RADIX_NLEN 2
#define M 1 // Number of window to check per iteration

mzd_t* semi_bc(mzd_t* G, uint64_t time_sec, uint64_t l, uint64_t c) {

    // Time mesuring stuff
    clock_t start = clock(), current;
    double elapsed = 0.0;

    // p is the Stern p parameter. (= number of LC to make for each rows subsets)
    uint64_t p = 2, iter = 0, mwin = 0, delta = 0, citer = 0, total = 0, done = 0, nb_total_cw = 0, offset = 0, lc_index = 0, idx = 0, delta_alt = 0;

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

    // nelem is k/2 choose p, number of elements in L1 and L2
    uint64_t nelem = ((N/4 * (N/4 - 1)) /2);

    // number of possible value on window
    uint64_t nb_keys = 1ULL << l; 


    // Holds the linear combinations made for each pass
    lc** lc_tab = (lc**)malloc(M * sizeof(lc*));
    lc** lc_tab_sorted = (lc**)malloc(M * sizeof(lc*));
    uint64_t** lc_offsets = (uint64_t**)malloc(M * sizeof(uint64_t*));
    CHECK_MALLOC(lc_tab);
    CHECK_MALLOC(lc_tab_sorted);
    CHECK_MALLOC(lc_offsets);

    for (mwin = 0; mwin < M; mwin++) {
        lc_tab[mwin] = (lc*)malloc(nelem * sizeof(lc));
        lc_tab_sorted[mwin] = (lc*)malloc(nelem * sizeof(lc));

        // lc_offsets[x] gives the offset at which we can find combination of value x in lc_tab
        lc_offsets[mwin] = (uint64_t*)malloc(nb_keys * sizeof(uint64_t));

        CHECK_MALLOC(lc_tab[mwin]);
        CHECK_MALLOC(lc_tab_sorted[mwin]);
    }

    // These arrays hold the size of their corresponding lc_tab
    uint64_t* lc_tab_size = (uint64_t*) malloc(M * sizeof(uint64_t));
    CHECK_MALLOC(lc_tab_size);

    // A generic holder for the current indexes
    uint64_t* lc_indexes = (uint64_t*)malloc(M * sizeof(uint64_t));
    CHECK_MALLOC(lc_indexes);

    // Generic alias on the sorted arrays
    lc** lc_tab_alias_sorted = (lc**)malloc(M * sizeof(lc*));

    // Precomputed mask for the window on which we want collision
    uint64_t sigmask = (1ULL << l) - 1;

    // Radixsort offsets array
    uint32_t *aux = (uint32_t*) malloc((1ULL << RADIX_WIDTH) * sizeof(uint32_t));

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
        memset(lc_indexes, 0, M * sizeof(uint64_t));

        // Creates all L1
        for (mwin = 0; mwin < M; mwin++) {
            for (comb1[0] = 0; comb1[0]  < 320 /* K/2 */; comb1[0]++) {
                uint64_t* row1 = (uint64_t*)Glw->rows[comb1[0]];

                for (comb1[1] = comb1[0] + 1; comb1[1] < 320 /* K/2 */; comb1[1]++) {
                    uint64_t* row2 = (uint64_t*)Glw->rows[comb1[1]];

                    // Compute the first l bits of the LC of rows 1 & 2 on the windows mwin
                    delta = (row1[mwin] ^ row2[mwin]) & sigmask;
                    lc_tab[mwin][lc_indexes[mwin]].index1 = comb1[0];
                    lc_tab[mwin][lc_indexes[mwin]].index2 = comb1[1];
                    lc_tab[mwin][lc_indexes[mwin]].delta = delta;
                    lc_indexes[mwin]++;
                }
            }
        }
    
        // Sort all L1
        for (mwin = 0; mwin < M; mwin++) {
            lc_tab_alias_sorted[mwin] = radixsort(lc_tab[mwin], lc_tab_sorted[mwin], nelem, aux);
        }

        // Compute offset tables for each L1
        for (mwin  = 0; mwin < M; mwin++) {
            memset(lc_offsets[mwin], 0, sizeof(uint32_t) * nb_keys);
            for (lc_index = 0; lc_index < nelem; lc_index++) {
                if (lc_offsets[mwin][lc_tab_alias_sorted[mwin][lc_index].delta] == 0) {
                    lc_offsets[mwin][lc_tab_alias_sorted[mwin][lc_index].delta] = lc_index;
                }
            }
        }
        
        // Creates all L2 on the fly
        for (comb2[0] = 320 /* K/2 */; comb2[0]  < 640 /* K */; comb2[0]++) {
            uint64_t* row1 = (uint64_t*)Glw->rows[comb2[0]];

            for (comb2[1] = comb2[0] + 1; comb2[1] < 640 /* K */; comb2[1]++) {
                uint64_t* row2 = (uint64_t*)Glw->rows[comb2[1]];

                for (mwin = 0; mwin < M; mwin++) {
                    // Compute the first l bits of the LC of rows 1 & 2
                    delta = (row1[mwin] ^ row2[mwin]) & sigmask;

                    for (idx = 0; idx < l + 1; idx++) {
                        delta_alt = delta;

                        // for idx in 0..l-1 we flip corresponding bit
                        // and for idx == l, we keep delta_alt = delta (exact match as in classical stern's algorithm)
                        if (idx < l) {
                            delta_alt ^= (1ULL << idx); // Flip bit number l
                        }
                    
                        offset = lc_offsets[mwin][delta_alt];

                        // Check that we have at least 1 result before building the combination of row1 + row2
                        if (lc_tab_alias_sorted[mwin][offset].delta == delta_alt) {

#if defined(AVX512_ENABLED)
                            __m512i linear_comb_high = _mm512_loadu_si512(row1);
                            __m128i linear_comb_low = _mm_loadu_si128(row1 + 64);

                            linear_comb_high = _mm512_xor_si512(linear_comb_high, _mm512_loadu_si512(row2));
                            linear_comb_low = _mm_xor_si128(linear_comb_low, _mm_loadu_si128(row2 + 64));
#else
                            memset(linear_comb, 0, 80); 
                            // Zero-out linear_comb and load row1 XOR row2 into it
                            mxor(linear_comb,(uint64_t*)linear_comb, 10);
                            mxor(linear_comb, (uint64_t*)Glw->rows[comb2[0]], 10);
                            mxor(linear_comb, (uint64_t*)Glw->rows[comb2[1]], 10);
#endif
                        
                            for (; offset < nelem && lc_tab_alias_sorted[mwin][offset].delta == delta_alt; offset++) {
                                comb1[0] = lc_tab_alias_sorted[mwin][offset].index1;
                                comb1[1] = lc_tab_alias_sorted[mwin][offset].index2;

#if defined(AVX512_ENABLED)
                                void* row3 = Glw->rows[comb1[0]];
                                void* row4 = Glw->rows[comb1[1]];

                                // Load the two new rows and add them to the LC of the two previous ones.
                                __m512i linear_comb_high_next = _mm512_xor_si512(linear_comb_high, _mm512_loadu_si512(row3));
                                __m128i linear_comb_low_next = _mm_xor_si128(linear_comb_low, _mm_loadu_si128(row3 + 64));

                                linear_comb_high_next = _mm512_xor_si512(linear_comb_high_next, _mm512_loadu_si512(row4));
                                linear_comb_low_next = _mm_xor_si128(linear_comb_low_next, _mm_loadu_si128(row4 + 64));

                                // Save the result of the LC of the 4 rows
                                _mm512_storeu_si512(linear_comb_next, linear_comb_high_next);
                                _mm_storeu_si128((__m128i*)(linear_comb_next + 8), linear_comb_low_next);
#else
                                memcpy(linear_comb_next, linear_comb, 80);
                                mxor(linear_comb_next, (uint64_t*)Glw->rows[comb1[0]], 10);
                                mxor(linear_comb_next, (uint64_t*)Glw->rows[comb1[1]], 10);
#endif

                                wt = popcnt64_unrolled(linear_comb_next , 10 /* K/64 */);
                                nb_total_cw++;

                                if (wt < min_wt) {

                                    // Save the new min weight and the indexes of th e linear combination to obtain it
                                    current = clock();
                                    elapsed = ((double)(current - start))/CLOCKS_PER_SEC;
                                    printf("niter=%lu, time=%.3f, wt=%ld\n", iter, elapsed, wt + 2*p);

                                    min_wt = wt;

                                    // Save the indexes of the LC
                                    memcpy(min_comb, comb1, 2 * sizeof(rci_t));
                                    memcpy(min_comb + 2, comb2, 2 * sizeof(rci_t));

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
    printf("# Average number of collisions / iter : %.3f\n", (double)nb_collisions/((double)iter * M));
    printf("# Average number of delta with at least 1 collision / nb delta : %.3f / %lu\n", (double)nb_collisions_delta/((double)iter * M), nelem);
    printf("# Average sorted list size : %.3f\n", (double)collisions_list_size/(2.0 * (double)iter * M));
#endif
    printf("# Total number of iterations done        %lu\n", iter);
    printf("# Iter/s                                 %.3f\n", (double)iter/(double)time_sec);
    printf("# Total number of analysed codewords     %lu\n", nb_total_cw);

    mzd_t* result =  stern_reconstruct_cw(min_comb, column_perms_copy, min_cw, p);

    mzd_free(Gtemp);
    mzd_free(Glw);

    free(min_cw);
    free(linear_comb);
    free(linear_comb_next);
    free(column_perms_copy);
    free(column_perms);

    for (mwin = 0; mwin < M; mwin++) {
        free(lc_tab[mwin]);
        free(lc_tab_sorted[mwin]);
    }

    free(lc_tab);
    free(lc_tab_sorted);


    free(lc_tab_size);
    free(lc_tab_alias_sorted);

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

lc* radixsort(lc* T, lc* Ts, int64_t Tlen, uint32_t* aux) {

    int i;
    lc* tmp;

    for (i = 0; i < RADIX_NLEN; i++) {
        tmp = T;
        T = denomsort_r(T, Ts, Tlen, RADIX_WIDTH, i*RADIX_WIDTH, aux);
        Ts = tmp;
    }

    return T;
}
