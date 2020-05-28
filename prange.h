#ifndef ISD_PRANGE_H
#define ISD_PRANGE_H

#if defined(__AVX512DQ__) && defined(__AVX512F__) && defined(__AVX512VL__)
#define AVX512_ENABLED
#endif

#include <m4ri.h>
#include <stdint.h>

// Input:
//  - G, a n/2 x n generator matrix
//  - niter, the number of iteration to make
//  Output:
//  - min_cw, the lowest codeword found
mzd_t* isd_prange_canteaut_chabaud(mzd_t* G, uint64_t niter);

#endif
