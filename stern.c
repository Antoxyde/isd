#include "stern.h"
#include "utils.h"
#include "xoshiro256starstar.h"
#include "libpopcnt.h"
#include <time.h>
#include <immintrin.h>



mzd_t* isd_stern_canteaut_chabaud_p2_sort(mzd_t* G, uint64_t niter, uint64_t sigma) {

    clock_t start = clock(), current;
    double elapsed = 0.0;

    uint64_t p = 2, iter = 0, nb_collision = 0;
    rci_t n = G->ncols, comb1[2], comb2[2], min_comb[4],lambda = 0, mu = 0, tmp = 0, i = 0, j = 0;
    rci_t k = n/2;

    int min_wt = 1000, wt = 0;
    void* row = NULL;
    (void)row; // otherwise gcc is :-(


#if defined(AVX512_ENABLED)
    uint64_t mask[10];
    memset(mask, 0, 80);
#endif

    uint64_t* word = NULL;
    uint64_t* linear_comb = (uint64_t*)malloc(sizeof(uint64_t) * 10);
    uint64_t* linear_comb_next = (uint64_t*)malloc(sizeof(uint64_t) * 10);
    uint64_t* min_cw = (uint64_t*)malloc(sizeof(uint64_t) * 10);
    rci_t* column_perms_copy =  (rci_t*) malloc(sizeof(rci_t) * n);
    rci_t* column_perms = (rci_t*) malloc(sizeof(rci_t) * n);
    CHECK_MALLOC(column_perms);
    CHECK_MALLOC(column_perms_copy);
    CHECK_MALLOC(linear_comb);
    CHECK_MALLOC(min_cw);

    for (i = 0; i < n; i++) column_perms[i] = i;

    mzd_t* Gtemp = mzd_copy(NULL, G);

    // Ensure that we work with a systematic generator matrix
    rref_to_systematic(Gtemp, column_perms);
    mzd_t* Glw = mzd_submatrix(NULL, Gtemp, 0, k, k, n);

    // Big array which contains all the linear combinations
    uint32_t nelem = ((n/4 * (n/4 - 1)) /2);
    lc* lc_tab = (lc*)malloc(nelem * sizeof(lc));
    lc* lc_tab_sorted = (lc*)malloc(nelem * sizeof(lc));

    lc* lc_tab_save = lc_tab;
    lc* lc_tab_sorted_save = lc_tab_sorted;

    uint32_t nb_keys = 1UL << sigma;
    uint32_t* lc_offsets = (uint32_t*)malloc(sizeof(uint32_t) * nb_keys);

    // Precomputed mask for the window on which we want collision
    uint64_t sigmask = (1 << sigma) - 1;

    //
    uint64_t radix_width = 9, radix_nlen = 2;
	uint32_t *aux = (uint32_t*) malloc((1 << radix_width) * sizeof(uint32_t));

    for (iter = 0; iter < niter; iter++) {

        // If we call radixsort with an odd nlen, lc_tab_sorted will point to lc_tab
        lc_tab_sorted = lc_tab_sorted_save;
        lc_tab = lc_tab_save;


        // Find lambda, mu s.t. Glw[lambda, mu] == 1
        lambda = xoshiro256starstar_random() % k;
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
        column_perms[mu + 640 /* k */] = tmp;

        // Clear the bit at the intersection of the lambda'th row and the mu'th column
        // so we don't have to rewrite a 1 in mu'th column everytime we add the lambda'th row
        mzd_write_bit(Glw, lambda, mu, 0);



#if defined(AVX512_ENABLED)
        void* row_lambda = Glw->rows[lambda];
        __m512i rlambda1 = _mm512_loadu_si512(row_lambda);
        __m128i rlambda2 = _mm_loadu_si128(row_lambda + 64 /* = 512 / (8 * sizeof(void)) */);

        // No easy instrinsic to set a single bit to 1 ?
        mask[big_mu] = (1ULL << small_mu);
        __m512i m1 = _mm512_loadu_si512(mask);
        __m128i m2 = _mm_loadu_si128(((void*)mask) + 64 /* = 512/(8 * sizeof(void)) */);
#endif

        // Add the lambda'th row to every other row that have a 1 in the mu'th column
        for (j = 0; j < k; j++) {
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

        // From here we have a new iset, we can start a Stern iteration
        uint64_t delta, lc_index = 0;

        // Gen all the LC from the first k/2 rows
        for (comb1[0] = 0; comb1[0]  < 320 /* n/4 */; comb1[0]++) {

            // Get the first sigma bits of the first row
            uint64_t row1 = ((uint64_t*)Glw->rows[comb1[0]])[0];

            for (comb1[1] = comb1[0] + 1; comb1[1] < n/4; comb1[1]++) {

                uint64_t row2 = ((uint64_t*)Glw->rows[comb1[1]])[0];

                // Compute the first sigma bits of the LC of rows 1 & 2
                delta = (row1 ^ row2) & sigmask;

                lc_tab[lc_index].index1 = comb1[0];
                lc_tab[lc_index].index2 = comb1[1];
                lc_tab[lc_index].delta = delta;
                lc_index++;
            }
        }

        //qsort(lc_tab, nelem, sizeof(lc), &compare_lc);
        lc_tab_sorted = radixsort(lc_tab, lc_tab_sorted, nelem, radix_width, radix_nlen, aux);

        lc_index = 0;

        memset(lc_offsets, 0, sizeof(uint32_t) * nb_keys);

        for (lc_index = 0; lc_index < nelem; lc_index++) {
            if (lc_offsets[lc_tab[lc_index].delta] == 0) {
                lc_offsets[lc_tab[lc_index].delta] = lc_index;
            }
        }


        for (comb2[0] = 320 /* n/4 */; comb2[0]  < 640 /* n/2 */; comb2[0]++) {

            uint64_t row1 = ((uint64_t*)Glw->rows[comb2[0]])[0];

            for (comb2[1] = comb2[0] + 1; comb2[1] < 640 /* n/2 */; comb2[1]++) {

                // Compute the "key" of the linear combination
                uint64_t row2 = ((uint64_t*)Glw->rows[comb2[1]])[0];

                // Compute the first sigma bits of the LC of rows 1 & 2
                delta = (row1 ^ row2) & sigmask;

                // And check if some elements from the previous set already had this key
                lc_index  = lc_offsets[delta];

#if defined(AVX512_ENABLED)
                void* row3 = Glw->rows[comb2[0]];
                void* row4 = Glw->rows[comb2[1]];

                __m512i linear_comb_high = _mm512_loadu_si512(row3);
                __m128i linear_comb_low = _mm_loadu_si128(row3 + 64);

                linear_comb_high = _mm512_xor_si512(linear_comb_high, _mm512_loadu_si512(row4));
                linear_comb_low = _mm_xor_si128(linear_comb_low, _mm_loadu_si128(row4 + 64));
#else

                mxor(linear_comb,(uint64_t*)linear_comb, 10);
                mxor(linear_comb, (uint64_t*)Glw->rows[comb2[0]], 10);
                mxor(linear_comb, (uint64_t*)Glw->rows[comb2[1]], 10);

#endif

                while (lc_index < nelem && lc_tab[lc_index].delta == delta) {

                    comb1[0] = lc_tab[lc_index].index1;
                    comb1[1] = lc_tab[lc_index].index2;

                    nb_collision++;
                    lc_index++;

#if defined(AVX512_ENABLED)
                    void* row1 = Glw->rows[comb1[0]];
                    void* row2 = Glw->rows[comb1[1]];

                    // Load the two new rows and add them to the LC of the two previous ones.
                    __m512i linear_comb_high_next = _mm512_xor_si512(linear_comb_high, _mm512_loadu_si512(row1));
                    __m128i linear_comb_low_next = _mm_xor_si128(linear_comb_low, _mm_loadu_si128(row1 + 64));

                    linear_comb_high_next = _mm512_xor_si512(linear_comb_high_next, _mm512_loadu_si512(row2));
                    linear_comb_low_next = _mm_xor_si128(linear_comb_low_next, _mm_loadu_si128(row2 + 64));

                    // Save the result of the LC of the 4 rows
                    _mm512_storeu_si512(linear_comb_next, linear_comb_high_next);
                    _mm_storeu_si128((__m128i*)(linear_comb_next + 8), linear_comb_low_next);
#else
                    memcpy(linear_comb_next, linear_comb, 80);
                    mxor(linear_comb_next, (uint64_t*)Glw->rows[comb1[0]], 10);
                    mxor(linear_comb_next, (uint64_t*)Glw->rows[comb1[1]], 10);

#endif
                    //printf("DBG Linear comb is : \n");
                    //printbin(linear_comb_next, 640);

                    wt = popcnt64_unrolled(linear_comb_next, 10);

                    if (wt < min_wt) {

                        // Save the new min weight and the indexes of th e linear combination to obtain it
                        current = clock();
                        elapsed = ((double)(current - start))/CLOCKS_PER_SEC;
                        printf("niter=%lu, time=%.3f, wt=%ld\n", iter, elapsed, wt + 2*p);

                        min_wt = wt;
                        // Save the indexes of the LC
                        memcpy(min_comb, comb1, p * sizeof(rci_t));
                        memcpy(min_comb + p, comb2, p * sizeof(rci_t));

                        memcpy(min_cw, linear_comb_next, 80);
                        memcpy(column_perms_copy, column_perms, n * sizeof(rci_t));
                    }
                }
            }
        }
    }

    // Reconstruct the codeword by concatenating the identity and unpermuting the columns
    mzd_t* ident = mzd_init(1, k);
    for (int i = 0; i < 2*p; i++)
        mzd_write_bit(ident, 0, min_comb[i], 1);

    mzd_t* min_cw_m = mzd_init(1, k);
    memcpy(mzd_first_row(min_cw_m), min_cw, 80);

    mzd_t* min_cw_full = mzd_concat(NULL, ident, min_cw_m);

    mzd_t* result = mzd_copy(NULL, min_cw_full);
    for (i = 0; i < n; i++) {
        if (i != column_perms_copy[i]) {
            mzd_write_bit(result, 0, column_perms_copy[i], mzd_read_bit(min_cw_full, 0, i));
        }
    }

    printf("nb_collision/iset : %.3f\n", (double)nb_collision/(double)niter);

    mzd_free(Gtemp);
    mzd_free(Glw);
    mzd_free(min_cw_full);
    mzd_free(min_cw_m);
    mzd_free(ident);

    free(min_cw);
    free(linear_comb);
    free(linear_comb_next);
    free(column_perms_copy);
    free(column_perms);
    free(lc_offsets);
    free(lc_tab);
	free(aux);

    return result;
}

int compare_lc(const void* a, const void* b) {
    return ((lc*)a)->delta > ((lc*)b)->delta;
}


lc* denomsort_r(lc* T, lc* Ts, int64_t Tlen, uint64_t width, uint64_t pos, uint32_t* Aux) {

	uint32_t mask, k;
	int64_t i;
	k = 1 << width;
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
        Aux[val]--;
		Ts[ Aux[val]] = T[i];
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
