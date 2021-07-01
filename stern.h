#ifndef ISD_STERN_H
#define ISD_STERN_H

#if defined(__AVX512DQ__) && defined(__AVX512F__) && defined(__AVX512VL__)
#define AVX512_ENABLED
#endif

#include "m4ri/m4ri.h"
#include <stdint.h>

// Represent a linear combination
typedef struct lc_ {
    uint16_t index1,index2;
    uint32_t delta;
} lc;


/* Input:
 *  - G, a n/2 x n generator matrix
 * - niter, the number of iteration to make
 * - sigma and p, the Stern parameters
 * - radix_width: the window size denomsort will use at each iteration
 * - radix_nlen: sigma/radix_width, number of iteration that radixsort will make.
 * - m, integer in {1,..,k/64}, the number of windows to force to zeroes
 * - c, intger >= 1, number of Gaussian elimination iterations to make to derive a new systematic generator matrix
 * - discard_nrows & discard_threshold : the function will first perform a popcnt on the first `discard_nwords` of a line, and if the hamming weight of that line is under `discard_threshold`, it will discard it
 *
 * Output:
 * - min_cw, the lowest codeword found
 */
 mzd_t* isd_stern_canteaut_chabaud_p2_sort(mzd_t* G, uint64_t niter, uint64_t sigma, uint64_t radix_width, uint64_t radix_nlen, uint64_t m, uint64_t c, uint64_t discard_threshold, uint64_t discard_nwords);

int compare_lc(const void* a, const void* b);
lc* denomsort_r(lc* T, lc* Ts, int64_t Tlen, uint64_t width, uint64_t pos, uint32_t* aux);
lc* radixsort(lc* T, lc* Ts, int64_t Tlen, uint64_t width, uint64_t nlen, uint32_t *aux);

#endif
