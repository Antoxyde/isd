#include "stern.h"
#include "utils.h"
#include "xoshiro256starstar.h"
#include "libpopcnt.h"
#include "combinations.h"
#include <time.h>
#include <immintrin.h>

mzd_t* stern(mzd_t* G, uint64_t time_sec) {

    // Time mesuring stuff
    clock_t start = clock(), current;
    double elapsed = 0.0;
    uint64_t iter = 0, delta = 0, cc_iter = 0, total = 0,  idx_win = 0, i = 0, j = 0;

    comb_t comb_struct;
    comb_diff_t comb_diff;

    uint16_t min_comb[P1+P2];

    int min_wt = K - 1, wt = 0;

    void* row = NULL;

#if defined(AVX512_ENABLED)
    uint64_t mask[10] /* K/64 */;
    memset(mask, 0, 80 /* K/8 */);
#endif /* AVX512_ENABLED */

#if defined(DEBUG)
    uint64_t nb_collisions = 0, nb_collisions_delta = 0, collisions_list_size = 0;
#endif /* DEBUG */

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
    
    // Init to identity permutation
    for (i = 0; i < N; i++) column_perms[i] = i;

    mzd_t* Gtemp = mzd_copy(NULL, G);

    // Ensure that we work with a systematic generator matrix
    rref_to_systematic(Gtemp, column_perms);
    mzd_t* Glw = mzd_submatrix(NULL, Gtemp, 0, K, K, N);

    // nelem is the max number of LCs we will make
    // so its an upper bound for the size of each array
    uint64_t nelem_p1 = binomial(320 /*K/2*/, P1);
    uint64_t nelem_p2 = binomial(320 /*K/2*/, P2);

    // Number of bits to store 1 for every possible different values for delta
    uint64_t nb_keys_bits = 1ULL << (L - 3);

#if defined(FILTER)
    // Used to track which combination are present
    uint64_t** collisions_first_pass = (uint64_t**)malloc( M * sizeof(uint64_t*));
    uint64_t** collisions_second_pass = (uint64_t**)malloc( M * sizeof(uint64_t*));
    CHECK_MALLOC(collisions_first_pass);
    CHECK_MALLOC(collisions_second_pass);
#endif /* FILTER */

    // Holds the linear combinations made for each pass
    lc_tab** L2s = (lc_tab**)malloc(M * sizeof(lc_tab*));
    lc_tab** L1s = (lc_tab**)malloc(M * sizeof(lc_tab*));
    lc_tab** L2s_sorted = (lc_tab**)malloc(M * sizeof(lc_tab*));
    lc_tab** L1s_sorted = (lc_tab**)malloc(M * sizeof(lc_tab*));

    CHECK_MALLOC(L2s);
    CHECK_MALLOC(L1s);
    CHECK_MALLOC(L2s_sorted);
    CHECK_MALLOC(L1s_sorted);

    for (idx_win = 0; idx_win < M; idx_win++) {

#if defined(FILTER)
        collisions_first_pass[idx_win] = (uint64_t*)malloc(nb_keys_bits);
        collisions_second_pass[idx_win] = (uint64_t*)malloc(nb_keys_bits);
        CHECK_MALLOC(collisions_first_pass[idx_win]);
        CHECK_MALLOC(collisions_second_pass[idx_win]);
#endif /* FILTER */

        L2s[idx_win] = malloc(sizeof(lc_tab));
        L1s[idx_win] = malloc(sizeof(lc_tab));
        L2s_sorted[idx_win] = malloc(sizeof(lc_tab));
        L1s_sorted[idx_win] = malloc(sizeof(lc_tab));

        CHECK_MALLOC(L2s[idx_win]);
        CHECK_MALLOC(L1s[idx_win]);
        CHECK_MALLOC(L2s_sorted[idx_win]);
        CHECK_MALLOC(L1s_sorted[idx_win]);
        
        L2s[idx_win]->p  = P2;
        L1s[idx_win]->p = P1;
        L2s_sorted[idx_win]->p  = P2;
        L1s_sorted[idx_win]->p = P1;

        L2s[idx_win]->current_size = 0;
        L1s[idx_win]->current_size = 0;
        L2s_sorted[idx_win]->current_size = 0;
        L1s_sorted[idx_win]->current_size = 0;

        L2s[idx_win]->lcs = LC_MALLOC(nelem_p2, P2); 
        L1s[idx_win]->lcs = LC_MALLOC(nelem_p1, P1);
        L2s_sorted[idx_win]->lcs = LC_MALLOC(nelem_p2, P2); 
        L1s_sorted[idx_win]->lcs = LC_MALLOC(nelem_p1, P1); 

        CHECK_MALLOC(L2s[idx_win]->lcs);
        CHECK_MALLOC(L1s[idx_win]->lcs);
        CHECK_MALLOC(L2s_sorted[idx_win]->lcs);
        CHECK_MALLOC(L1s_sorted[idx_win]->lcs);
    }

    // A generic holder for the current indexes
    uint64_t* lc_indexes = (uint64_t*)malloc(M * sizeof(uint64_t));
    CHECK_MALLOC(lc_indexes);

    // Generic alias on the sorted arrays

    lc_tab** L2s_alias_sorted = (lc_tab**)malloc(M * sizeof(lc_tab*));
    lc_tab** L1s_alias_sorted = (lc_tab**)malloc(M * sizeof(lc_tab*));

    // Precomputed mask for the window on which we want collision
    uint64_t sigmask = (1ULL << L) - 1;

    // Radixsort offsets array
    uint32_t *aux = (uint32_t*) malloc((1ULL << RADIX_WIDTH) * sizeof(uint32_t));
    CHECK_MALLOC(aux);
    
    // Combinations tabs for enumeration of pX choose k/2 elements
    uint16_t* combinations_p1 = (uint16_t*)malloc((P1 + 1) * sizeof(uint16_t));
    uint16_t* combinations_p2 = (uint16_t*)malloc((P2 + 1) * sizeof(uint16_t));
    CHECK_MALLOC(combinations_p1);
    CHECK_MALLOC(combinations_p2);


    while (1) {

        /* Start of the Canteaut-Chabaud stuff, to derive a new information set. */
        for (cc_iter = 0; cc_iter < CC; cc_iter++) {
            canteaut_chabaud(Glw, column_perms);
        }
        /* End of the Canteaut-Chabaud stuff. We now have a proper Iset to work with. */
        
        // Reset previous iteration stuff
        for (idx_win = 0; idx_win < M; idx_win++) {

#if defined(FILTER)
            memset(collisions_first_pass[idx_win], 0, nb_keys_bits);
            memset(collisions_second_pass[idx_win], 0, nb_keys_bits);
#endif /* FILTER */

            L2s[idx_win]->current_size = 0;
            L1s[idx_win]->current_size = 0;
            L2s_sorted[idx_win]->current_size = 0;
            L1s_sorted[idx_win]->current_size = 0;
        }

#if defined(FILTER)
        // 1st pass
        memset(combinations_p1, 0, (P1 + 1) * sizeof(uint16_t));
        init_combination(&comb_struct,combinations_p1, P1, 320 /*K/2*/);

        for (i = 0; i < nelem_p1; i++) {
            for (idx_win = 0; idx_win < M; idx_win++) {
                delta = 0;
                for (j = 0; j < P1; j++) {
                    delta ^= Glw->rows[comb_struct.combination[j]][idx_win];
                }

                delta &= sigmask;

                STERN_SET_ONE(collisions_first_pass[idx_win], delta);
            }
            
            next_combination(&comb_struct, &comb_diff);
        }
#endif /* FILTER */
        
        // 2nd pass
        memset(combinations_p2, 0, (P2 + 1) * sizeof(uint16_t));
        init_combination(&comb_struct,combinations_p2, P2, 320 /*K/2*/);
        for (i = 0; i < nelem_p2; i++) {
            for (idx_win = 0; idx_win < M; idx_win++) {

                delta = 0;
                for (j = 0; j < P2; j++) {
                    delta ^= Glw->rows[comb_struct.combination[j] + 320][idx_win];
                }

                delta &= sigmask;

#if defined(FILTER)
                if (STERN_GET(collisions_first_pass[idx_win], delta)) {
#endif /* FILTER */

                    lc* elem = LC_TAB_GET(L2s[idx_win], L2s[idx_win]->current_size);
                    memcpy(elem->indexes, comb_struct.combination, P2 * sizeof(uint16_t));

                    elem->delta = delta;
                    L2s[idx_win]->current_size++;
#if defined(FILTER)
                    STERN_SET_ONE(collisions_second_pass[idx_win], delta);
                }
#endif /* FILTER */
            }
            next_combination(&comb_struct, &comb_diff);
        }
        
        // 3rd pass
        memset(combinations_p1, 0, (P1 + 1) * sizeof(uint16_t));
        init_combination(&comb_struct,combinations_p1, P1, 320 /*K/2*/);
        for (i = 0; i < nelem_p1; i++) {
            for (idx_win = 0; idx_win < M; idx_win++) {

                delta = 0;
                for (j = 0; j < P1; j++) {
                    delta ^= Glw->rows[comb_struct.combination[j]][idx_win];
                }

                delta &= sigmask;
#if defined(FILTER)
                if (STERN_GET(collisions_second_pass[idx_win], delta)) {
#endif /* FILTER */
                    lc* elem = LC_TAB_GET(L1s[idx_win], L1s[idx_win]->current_size);
                    memcpy(elem->indexes, comb_struct.combination, P1 * sizeof(uint16_t));
                    elem->delta = delta;
                    L1s[idx_win]->current_size++;
#if defined(FILTER)
                }
#endif /* FILTER */

            }
            next_combination(&comb_struct, &comb_diff);
        }

        /* From here, L1s[i][:L1s_size[i]] and L2s[i][:L2s_size[i]]
        *  contains the LCs that WILL collide.
        * So we sort em all. */


        for (idx_win = 0; idx_win < M; idx_win++) {
            
            L2s_sorted[idx_win]->current_size = L2s[idx_win]->current_size;
            L1s_sorted[idx_win]->current_size = L1s[idx_win]->current_size;

            // We have to use aliases to avoid mixing up the LXs_sorted and LXs pointers (radixsort returns one or the other, depending on the radix_nlen parameter.)
            L2s_alias_sorted[idx_win] = radixsort(L2s[idx_win], L2s_sorted[idx_win], L2s[idx_win]->current_size, aux);
            L1s_alias_sorted[idx_win] = radixsort(L1s[idx_win], L1s_sorted[idx_win], L1s[idx_win]->current_size, aux);
        }
        
        // Now that everything is sorted,  we can actually match all the pairs to get the collisions.
        for (idx_win = 0; idx_win < M; idx_win++) {

            int index_second = 0, index_third = 0, save_index_third = 0;

            for (index_second = 0; index_second < L2s_alias_sorted[idx_win]->current_size; index_second++) {
                
                if (index_second > 0 && LC_TAB_GET(L2s_alias_sorted[idx_win],index_second - 1)->delta != LC_TAB_GET(L2s_alias_sorted[idx_win],index_second)->delta) {
                    save_index_third = index_third;
                }

                // Xor the first P2 rows
#if defined(AVX512_ENABLED)
                row = (void*)Glw->rows[LC_TAB_GET(L2s_alias_sorted[idx_win], index_second)->indexes[0] + 320];

                __m512i linear_comb_high = _mm512_loadu_si512(row);
                __m128i linear_comb_low = _mm_loadu_si128(row + 64 /* 512/8 */ );
                for (i = 1; i < P2; i++) {
                    row = (void*)Glw->rows[LC_TAB_GET(L2s_alias_sorted[idx_win], index_second)->indexes[i] + 320];
                    linear_comb_high = _mm512_xor_si512(linear_comb_high, _mm512_loadu_si512(row));
                    linear_comb_low = _mm_xor_si128(linear_comb_low, _mm_loadu_si128(row + 64 /* 512/8 */));
                }
#else
                memset(linear_comb, 0, 80 /* K/8 */);
                for (i = 0; i < P2; i++) {
                    mxor(linear_comb, (uint64_t*)Glw->rows[LC_TAB_GET(L2s_alias_sorted[idx_win], index_second)->indexes[i] + 320], 10 /* K/64 */);
                }
#endif /* AVX512_ENABLED */
                for (index_third = save_index_third; index_third < L1s_alias_sorted[idx_win]->current_size && LC_TAB_GET(L2s_alias_sorted[idx_win], index_second)->delta >= LC_TAB_GET(L1s_alias_sorted[idx_win], index_third)->delta; index_third++) {

                    if (LC_TAB_GET(L2s_alias_sorted[idx_win], index_second)->delta > LC_TAB_GET(L1s_alias_sorted[idx_win], index_third)->delta) continue;

#if defined(DEBUG)
                    nb_collisions++;
#endif /* DEBUG */

                    // Xor the P1 other rows and xor it to the first result
#if defined(AVX512_ENABLED)
                    row = (void*)Glw->rows[LC_TAB_GET(L1s_alias_sorted[idx_win],index_third)->indexes[0]];

                    // Load the two new rows and add them to the LC of the two previous ones.
                    __m512i linear_comb_high_next = _mm512_xor_si512(linear_comb_high, _mm512_loadu_si512(row));
                    __m128i linear_comb_low_next = _mm_xor_si128(linear_comb_low, _mm_loadu_si128(row + 64 /* 512/(8 * sizeof(void)) */));

                    for (i = 1; i < P1; i++) {
                        row = (void*)Glw->rows[LC_TAB_GET(L1s_alias_sorted[idx_win], index_second)->indexes[i]];
                        linear_comb_high_next = _mm512_xor_si512(linear_comb_high_next, _mm512_loadu_si512(row));
                        linear_comb_low_next = _mm_xor_si128(linear_comb_low_next, _mm_loadu_si128(row + 64 /* 512/(8 * sizeof(void)) */));
                    }

                    // Save the result of the LC of the P1+P2 rows
                    _mm512_storeu_si512(linear_comb_next, linear_comb_high_next);
                    _mm_storeu_si128((__m128i*)(linear_comb_next + 8 /* 512/(8 * sizeof(uint64_t) */), linear_comb_low_next);

#else
                    memcpy(linear_comb_next, linear_comb, 80 /* K/8 */ );
                    for (i = 0; i < P1; i++) {
                        mxor(linear_comb_next, (uint64_t*)Glw->rows[LC_TAB_GET(L1s_alias_sorted[idx_win],index_third)->indexes[i]], 10 /* K/64 */);
                    }
#endif /* AVX512_ENABLED */

                    wt = popcnt64_unrolled(linear_comb_next, 10 /* K/64 */);
                    total++;
                    if (wt < min_wt) {

                        // Save the new min weight and the indexes of th e linear combination to obtain it
                        current = clock();
                        elapsed = ((double)(current - start))/CLOCKS_PER_SEC;
                        printf("niter=%lu, time=%.3f, wt=%d\n", iter, elapsed, wt + P1 + P2);

                        min_wt = wt;
                        // Save the indexes of the LC 
                        memcpy(min_comb, LC_TAB_GET(L2s_alias_sorted[idx_win], index_second)->indexes, P2 * sizeof(uint16_t));
                        for (i = 0; i < P2; i++) {
                            min_comb[i] += 320;
                        }

                        memcpy(min_comb + P2, LC_TAB_GET(L1s_alias_sorted[idx_win], index_third)->indexes, P1 * sizeof(uint16_t));

                        memcpy(min_cw, linear_comb_next, 80 /* K/8 */);
                        memcpy(column_perms_copy, column_perms, N * sizeof(rci_t));
                        
  
                        mzd_t* cw = stern_reconstruct_cw(min_comb, column_perms_copy, min_cw, P1 + P2);
                        print_cw(cw);
                        mzd_free(cw);

                        fflush(stdout);
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

#if defined(DEBUG)
    printf("# Average number of collisions / iter : %.3f\n", (double)nb_collisions/((double)iter * M));
#endif /* DEBUG */
    printf("# Total number of iterations done : %lu\n", iter);
    printf("# Iter/s : %.3f\n", (double)iter/(double)time_sec);
    printf("# Number of words analysed : %lu\n", total);

    mzd_t* result =  stern_reconstruct_cw(min_comb, column_perms_copy, min_cw, P1 + P2);

    mzd_free(Gtemp);
    mzd_free(Glw);

    free(min_cw);
    free(linear_comb);
    free(linear_comb_next);
    free(column_perms_copy);
    free(column_perms);

    for (idx_win = 0; idx_win < M; idx_win++) {


        free(L2s[idx_win]->lcs);
        free(L1s[idx_win]->lcs);
        free(L2s_sorted[idx_win]->lcs);
        free(L1s_sorted[idx_win]->lcs);
    }

    free(L2s);
    free(L1s);
    free(L2s_sorted);
    free(L1s_sorted);

    free(L2s_alias_sorted);
    free(L1s_alias_sorted);

    free(lc_indexes);
    free(aux);

#if defined(FILTER)
    for (idx_win = 0; idx_win < M; idx_win++) {
        free(collisions_first_pass[idx_win]);
        free(collisions_second_pass[idx_win]);
    }

    free(collisions_first_pass);
    free(collisions_second_pass);
#endif /* FILTER */


    return result;
}


lc_tab* denomsort_r(lc_tab* T, lc_tab* Ts, int64_t Tlen, uint64_t width, uint64_t pos, uint32_t* Aux) {

    uint32_t mask, k;
    int64_t i;
    k = 1UL << width;
    mask = k - 1;

    memset(Aux, 0, k * sizeof(uint32_t));

    for (i = 0; i < Tlen; i++) {
        Aux[ ( LC_TAB_GET(T, i)->delta >> pos) & mask]++;
    }

    for (i = 1; i < k; i++) {
        Aux[i] += Aux[i - 1];
    }

    for (i = Tlen - 1; i >= 0; i--) {
        uint32_t val = (LC_TAB_GET(T, i)->delta >> pos) & mask;
        Aux[val]--; 
        memcpy(LC_TAB_GET(Ts, Aux[val]), LC_TAB_GET(T, i), sizeof(lc) + T->p * sizeof(uint16_t));
    }

    return Ts;
}

lc_tab* radixsort(lc_tab* T, lc_tab* Ts, int64_t Tlen, uint32_t* aux) {

    int i;
    lc_tab* tmp;

    for (i = 0; i < RADIX_LEN; i++) {
        tmp = T;
        T = denomsort_r(T, Ts, Tlen, RADIX_WIDTH, i*RADIX_WIDTH, aux);
        Ts = tmp;
    }

    return T;
}
