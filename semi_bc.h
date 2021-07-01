#ifndef ISD_SEMIBC_H
#define ISD_SEMIBC_H

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
 * - l , the Stern's algorithm parameter for window size
 * - c, intger >= 1, number of Gaussian elimination iterations to make to derive a new systematic generator matrix
 *
 * Output:
 * - min_cw, the lowest codeword found
 *
 *   Stern's algorithm parameter p is fixed to 2
 */
 mzd_t* semi_bc(mzd_t* G, uint64_t niter, uint64_t l,uint64_t c);

int compare_lc(const void* a, const void* b);
lc* denomsort_r(lc* T, lc* Ts, int64_t Tlen, uint64_t width, uint64_t pos, uint32_t* aux);
lc* radixsort(lc* T, lc* Ts, int64_t Tlen, uint32_t *aux);

#endif
