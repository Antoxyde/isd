#ifndef ISD_ISD_H
#define ISD_ISD_H

#include <m4ri.h>

// Input:
//  - G, a n/2 x n generator matrix
//  - niter, the number of iteration to make
//  Output:
//  - min_cw, the lowest codeword found
mzd_t* isd_prange(mzd_t* G, int niter);

// Input:
//  - G, a n/2 x n generator matrix
//  - niter, the number of iteration to make
//  Output:
//  - min_cw, the lowest codeword found
mzd_t* isd_prange_canteaut(mzd_t* G, int niter);

// Input:
//  - G, a n/2 x n generator matrix
//  - niter, the number of iteration to make
//  Output:
//  - min_cw, the lowest codeword found
mzd_t* isd_prange_canteaut_test(mzd_t* G, int niter);


#endif
