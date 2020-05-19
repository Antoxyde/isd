#ifndef ISD_ISD_H
#define ISD_ISD_H

#include <m4ri.h>

// Input:
//  - G, an n/2 x n generator matrix
//  - Gis, a n/2 x n/2 matrix, in which we will store the random iset
//  - indices, a set of values {0..n} to be shuffled
void get_random_iset(const mzd_t* Gt, mzd_t* Gis, mzd_t* Gist, rci_t* indices);


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
mzd_t* isd_lee_brickell(mzd_t* G, int niter);

#endif
