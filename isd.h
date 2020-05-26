#ifndef ISD_ISD_H
#define ISD_ISD_H

#if defined(__AVX512DQ__) && defined(__AVX512F__) && defined(__AVX512VL__)
#define HAVE_AVX512
#endif

#include <m4ri.h>

// Input:
//  - G, a n/2 x n generator matrix
//  - niter, the number of iteration to make
//  Output:
//  - min_cw, the lowest codeword found
mzd_t* isd_prange_canteaut(mzd_t* G, int niter);


#endif
