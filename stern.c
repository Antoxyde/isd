#include "stern.h"
#include "utils.h"
#include "xoshiro256starstar.h"
#include "libpopcnt.h"
#include <time.h>
#include <immintrin.h>

mzd_t* stern(mzd_t* G, uint64_t time_sec) {

    // Time mesuring stuff
    clock_t start = clock(), current;
    double elapsed = 0.0;
    uint64_t iter = 0, delta = 0, cc_iter = 0, total = 0, done = 0, idx_win = 0;

    rci_t comb1[P1], comb2[P2], min_comb[P1+P2], i = 0;

    int min_wt = K - 1, wt = 0;

#if defined(AVX512_ENABLED)
    uint64_t mask[10] /* K/64 */;
    memset(mask, 0, 80 /* K/8 */);
#endif

#if defined(DEBUG)
    uint64_t nb_collisions = 0, nb_collisions_delta = 0, collisions_list_size = 0;
#endif

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
    uint64_t nelem = ((N/4 * (N/4 - 1)) /2); // TODO

    // Number of bits to store 1 for every possible different values for delta
    uint64_t nb_keys_bits = 1ULL << (L - 3);

    // Used to track which combination are present
    uint64_t** collisions_first_pass = (uint64_t**)malloc( M * sizeof(uint64_t*));
    uint64_t** collisions_second_pass = (uint64_t**)malloc( M * sizeof(uint64_t*));

    // Holds the linear combinations made for each pass
    lc_tab** L2s = (lc_tab**)malloc(M * sizeof(lc_tab*));
    lc_tab** L3s = (lc_tab**)malloc(M * sizeof(lc_tab*));
    lc_tab** L2s_sorted = (lc_tab**)malloc(M * sizeof(lc_tab*));
    lc_tab** L3s_sorted = (lc_tab**)malloc(M * sizeof(lc_tab*));

    CHECK_MALLOC(collisions_first_pass);
    CHECK_MALLOC(collisions_second_pass);
    CHECK_MALLOC(L2s);
    CHECK_MALLOC(L3s);
    CHECK_MALLOC(L2s_sorted);
    CHECK_MALLOC(L3s_sorted);

    for (idx_win = 0; idx_win < M; idx_win++) {
        collisions_first_pass[idx_win] = (uint64_t*)malloc(nb_keys_bits);
        collisions_second_pass[idx_win] = (uint64_t*)malloc(nb_keys_bits);

        L2s[idx_win] = malloc(sizeof(lc_tab));
        L3s[idx_win] = malloc(sizeof(lc_tab));
        L2s_sorted[idx_win] = malloc(sizeof(lc_tab));
        L3s_sorted[idx_win] = malloc(sizeof(lc_tab));

        CHECK_MALLOC(L2s[idx_win]);
        CHECK_MALLOC(L3s[idx_win]);
        CHECK_MALLOC(L2s_sorted[idx_win]);
        CHECK_MALLOC(L3s_sorted[idx_win]);
        
        L2s[idx_win]->p  = P2;
        L3s[idx_win]->p = P1;
        L2s_sorted[idx_win]->p  = P2;
        L3s_sorted[idx_win]->p = P1;

        L2s[idx_win]->lcs = LC_MALLOC(nelem, P2); 
        L3s[idx_win]->lcs = LC_MALLOC(nelem, P1);
        L2s_sorted[idx_win]->lcs = LC_MALLOC(nelem, P2); 
        L3s_sorted[idx_win]->lcs = LC_MALLOC(nelem, P1); 

        CHECK_MALLOC(collisions_first_pass[idx_win]);
        CHECK_MALLOC(collisions_second_pass[idx_win]);
        CHECK_MALLOC(L2s[idx_win]->lcs);
        CHECK_MALLOC(L3s[idx_win]->lcs);
        CHECK_MALLOC(L2s_sorted[idx_win]->lcs);
        CHECK_MALLOC(L3s_sorted[idx_win]->lcs);
    }

    // Holds the size of their corresponding lc_tab
    uint64_t* L2s_size = (uint64_t*) malloc(M * sizeof(uint64_t));
    uint64_t* L3s_size = (uint64_t*) malloc(M * sizeof(uint64_t));
    CHECK_MALLOC(L2s_size);
    CHECK_MALLOC(L3s_size);

    // A generic holder for the current indexes
    uint64_t* lc_indexes = (uint64_t*)malloc(M * sizeof(uint64_t));
    CHECK_MALLOC(lc_indexes);

    // Generic alias on the sorted arrays

    lc_tab** L2s_alias_sorted = (lc_tab**)malloc(M * sizeof(lc_tab*));
    lc_tab** L3s_alias_sorted = (lc_tab**)malloc(M * sizeof(lc_tab*));

    // Precomputed mask for the window on which we want collision
    uint64_t sigmask = (1ULL << L) - 1;

    // Radixsort offsets array
    uint32_t *aux = (uint32_t*) malloc((1ULL << RADIX_WIDTH) * sizeof(uint32_t));

    while (1) {

        /* Start of the Canteaut-Chabaud stuff, to derive a new information set. */
        for (cc_iter = 0; cc_iter < CC; cc_iter++) {
            canteaut_chabaud(Glw, column_perms);
        }
        /* End of the Canteaut-Chabaud stuff. We now have a proper Iset to work with. */

        // Reset all the stuff we will need in the iteration
        memset(lc_indexes, 0, M * sizeof(uint64_t));

        for (idx_win = 0; idx_win < M; idx_win++) {
            memset(collisions_first_pass[idx_win], 0, nb_keys_bits);
            memset(collisions_second_pass[idx_win], 0, nb_keys_bits);
        }
        

        // TODO variable P1
        // 1st pass, gen all the LC from the first k/2 rows
        for (comb1[0] = 0; comb1[0]  < 320 /* K/2 */; comb1[0]++) {
            uint64_t* row1 = (uint64_t*)Glw->rows[comb1[0]];

            for (comb1[1] = comb1[0] + 1; comb1[1] < 320 /* K/2 */; comb1[1]++) {
                uint64_t* row2 = (uint64_t*)Glw->rows[comb1[1]];

                for (idx_win = 0; idx_win < M; idx_win++) {
                    // Compute the first sigma bits of the LC of rows 1 & 2 on the windows idx_win
                    delta = (row1[idx_win] ^ row2[idx_win]) & sigmask;
                    STERN_SET_ONE(collisions_first_pass[idx_win], delta);
                }

            }
        }
    
        // TODO variable P2
        // 2nd pass, gen all the LC from the k/2 to k rows but store
        // only the ones that will collides with at least 1 element from lc_tab_first
        for (comb2[0] = 320 /* K/2 */; comb2[0]  < 640 /* K */; comb2[0]++) {

            uint64_t* row1 = (uint64_t*)Glw->rows[comb2[0]];

            for (comb2[1] = comb2[0] + 1; comb2[1] < 640 /* K */; comb2[1]++) {

                uint64_t* row2 = (uint64_t*)Glw->rows[comb2[1]];

                // The `idx_win`st window is the first sigma bits of the `idx_win`th word of the row.
                for (idx_win = 0; idx_win < M; idx_win++) {

                    // Compute the first sigma bits of the LC of rows 1 & 2
                    delta = (row1[idx_win] ^ row2[idx_win]) & sigmask;

                    if (STERN_GET(collisions_first_pass[idx_win], delta)) {
                        lc* elem = LC_TAB_GET(L2s[idx_win], lc_indexes[idx_win]);
                        elem->indexes[0] = comb2[0];
                        elem->indexes[1] = comb2[1];
                        elem->delta = delta;
                        STERN_SET_ONE(collisions_second_pass[idx_win], delta);
                        lc_indexes[idx_win]++;
                    }
                }
            }
        }

        // Save the size of each lc_tab_first elements and reset lc_indexes
        memcpy(L2s_size, lc_indexes, M * sizeof(uint64_t));
        memset(lc_indexes, 0, M * sizeof(uint64_t));

        // 3rd pass, copy all the element of first tab that will actually collide in a new array
        for (comb1[0] = 0; comb1[0]  < 320 /* K/2 */; comb1[0]++) {
            uint64_t* row1 = (uint64_t*)Glw->rows[comb1[0]];

            for (comb1[1] = comb1[0] + 1; comb1[1] < 320 /* K/2 */; comb1[1]++) {
                uint64_t* row2 = (uint64_t*)Glw->rows[comb1[1]];

                for (idx_win = 0; idx_win < M; idx_win++) {

                    delta = (row1[idx_win] ^ row2[idx_win]) & sigmask;

                    if (STERN_GET(collisions_second_pass[idx_win], delta)) {

                        lc* elem = LC_TAB_GET(L3s[idx_win], lc_indexes[idx_win]);
                        elem->indexes[0] = comb1[0];
                        elem->indexes[1] = comb1[1];
                        elem->delta = delta;
                        lc_indexes[idx_win]++;
                    }
                }
            }
        }


        // Save the size of each lc_tab_first elements and reset lc_indexes
        memcpy(L3s_size, lc_indexes, M * sizeof(uint64_t));

        /* From here, L3s[i][:L3s_size[i]] and L2s[i][:L2s_size[i]]
        *  contains the LCs that WILL collide.
        * So we sort em all. */

        for (idx_win = 0; idx_win < M; idx_win++) {

#if defined(DEBUG)
            collisions_list_size += L2s_size[idx_win] + L3s_size[idx_win];
#endif

            // We have to use aliases to avoid mixing up the lc_tab_X and lc_tab_X_sorted pointers (radixsort returns one of them, depending on the radix_nlen parameter.)
            L2s_alias_sorted[idx_win] = radixsort(L2s[idx_win], L2s_sorted[idx_win], L2s_size[idx_win], aux);
            L3s_alias_sorted[idx_win] = radixsort(L3s[idx_win], L3s_sorted[idx_win], L3s_size[idx_win], aux);
        }
        
        /*
        printf("Sorted : \n");
        for (int mwin = 0; mwin < M; mwin++) {
            for (i = 0; i < 15; i++) {
                lc* e = LC_TAB_GET(L2s_alias_sorted[mwin], i);
                lc* e2 = LC_TAB_GET(L2s_alias_sorted[mwin], i+1);
                printf("%d : %d [%d,%d]\n", i, e->delta, e->indexes[0], e->indexes[1]);
                printf("%d +1: %d [%d,%d]\n", i, e2->delta, e2->indexes[0], e2->indexes[1]);
            }
        }
        */


        // Now that everything is sorted,  we can actually match all the pairs to get the collisions.
        for (idx_win = 0; idx_win < M; idx_win++) {

            int index_second = 0, index_third = 0, save_index_third = 0;

            for (index_second = 0; index_second < L2s_size[idx_win]; index_second++) {
                
                if (index_second > 0 && LC_TAB_GET(L2s_alias_sorted[idx_win],index_second - 1)->delta != LC_TAB_GET(L2s_alias_sorted[idx_win],index_second)->delta) {
                    save_index_third = index_third;

#if defined(DEBUG)
                    nb_collisions_delta++;
#endif
                }

                // Xor the first 2 rows
#if defined(AVX512_ENABLED)
                void* row1 = (void*)Glw->rows[LC_TAB_GET(L2s_alias_sorted[idx_win], index_second)->indexes[0]];
                void* row2 = (void*)Glw->rows[LC_TAB_GET(L2s_alias_sorted[idx_win], index_second)->indexes[1]];

                __m512i linear_comb_high = _mm512_loadu_si512(row1);
                __m128i linear_comb_low = _mm_loadu_si128(row1 + 64 /* 512/8 */ );

                linear_comb_high = _mm512_xor_si512(linear_comb_high, _mm512_loadu_si512(row2));
                linear_comb_low = _mm_xor_si128(linear_comb_low, _mm_loadu_si128(row2 + 64 /* 512/8 */));
#else
                memset(linear_comb, 0, 80 /* K/8 */);
                mxor(linear_comb, (uint64_t*)Glw->rows[LC_TAB_GET(L2s_alias_sorted[idx_win], index_second)->indexes[0]], 10 /* K/64 */);
                mxor(linear_comb, (uint64_t*)Glw->rows[LC_TAB_GET(L2s_alias_sorted[idx_win],index_second)->indexes[1]], 10 /* K/64 */);
#endif

                for (index_third = save_index_third; index_third < L3s_size[idx_win] && LC_TAB_GET(L2s_alias_sorted[idx_win], index_second)->delta == LC_TAB_GET(L3s_alias_sorted[idx_win], index_third)->delta; index_third++) {

#if defined(DEBUG)
                    nb_collisions++;
#endif

                    // Xor the two other rows and xor it to the first result
#if defined(AVX512_ENABLED)
                    void* row3 = (void*)Glw->rows[LC_TAB_GET(L3s_alias_sorted[idx_win],index_third)->indexes[0]];
                    void* row4 = (void*)Glw->rows[LC_TAB_GET(L3s_alias_sorted[idx_win],index_third)->indexes[1]];

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
                    mxor(linear_comb_next, (uint64_t*)Glw->rows[LC_TAB_GET(L3s_alias_sorted[idx_win],index_third)->indexes[0]], 10 /* K/64 */);
                    mxor(linear_comb_next, (uint64_t*)Glw->rows[LC_TAB_GET(L3s_alias_sorted[idx_win],index_third)->indexes[1]], 10 /* K/64 */);
#endif

                    //printf("DBG idx_win = %lu, linear comb is : \n", idx_win);
                    //printbin(linear_comb_next, 640);

                    wt = popcnt64_unrolled(linear_comb_next, 10 /* K/64 */);
                    total++;
                    if (wt < min_wt) {

                        // Save the new min weight and the indexes of th e linear combination to obtain it
                        current = clock();
                        elapsed = ((double)(current - start))/CLOCKS_PER_SEC;
                        printf("niter=%lu, time=%.3f, wt=%d\n", iter, elapsed, wt + P1 + P2);

                        min_wt = wt;

                        // Save the indexes of the LC TODO modularisé P1 P2
                        min_comb[0] = LC_TAB_GET(L2s_alias_sorted[idx_win], index_second)->indexes[0];
                        min_comb[1] = LC_TAB_GET(L2s_alias_sorted[idx_win], index_second)->indexes[1];
                        min_comb[2] = LC_TAB_GET(L3s_alias_sorted[idx_win], index_third)->indexes[0];
                        min_comb[3] = LC_TAB_GET(L3s_alias_sorted[idx_win], index_third)->indexes[1];

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

#ifdef DEBUG
    printf("# Average number of collisions / iter : %.3f\n", (double)nb_collisions/((double)iter * M));
    printf("# Average number of delta with at least 1 collision / nb delta : %.3f / %lu\n", (double)nb_collisions_delta/((double)iter * M), nelem);
    printf("# Average sorted list size : %.3f\n", (double)collisions_list_size/(2.0 * (double)iter * M));
#endif
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
        free(collisions_first_pass[idx_win]);
        free(collisions_second_pass[idx_win]);
        free(L2s[idx_win]->lcs);
        free(L3s[idx_win]->lcs);
        free(L2s_sorted[idx_win]->lcs);
        free(L3s_sorted[idx_win]->lcs);
    }

    free(L2s);
    free(L3s);
    free(L2s_sorted);
    free(L3s_sorted);

    free(collisions_first_pass);
    free(collisions_second_pass);

    free(L2s_size);
    free(L3s_size);
    free(L2s_alias_sorted);
    free(L3s_alias_sorted);

    free(lc_indexes);

    free(aux);

    return result;
}


lc_tab* denomsort_r(lc_tab* T, lc_tab* Ts, int64_t Tlen, uint64_t width, uint64_t pos, uint32_t* Aux) {

    uint32_t mask, k;
    int64_t i;
    k = 1UL << width;
    mask = k - 1;

    memset(Aux, 0, k * sizeof(uint32_t));

    for (i = 0; i < Tlen; i++) {
        lc* e = LC_TAB_GET(T, i);
        Aux[ (e->delta >> pos) & mask]++;
    }

    for (i = 1; i < k; i++) {
        Aux[i] += Aux[i - 1];
    }

    for (i = Tlen - 1; i >= 0; i--) {
        lc* e = LC_TAB_GET(T, i);
        uint32_t val = (e->delta >> pos) & mask;
        Aux[val]--; 
        lc* out = LC_TAB_GET(Ts, Aux[val]);
        lc* in = LC_TAB_GET(T, i);
        *out = *in;
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
