#ifndef ISD_ISD_H
#define ISD_ISD_H

#include <m4ri.h>

// Input:
//  - G, a n/2 x n generator matrix
//  - niter, the number of iteration to make
//
//  Output:
//  - min_cw, the lowest codeword found
mzd_t* isd_prange(mzd_t* G, int niter);

// Input:
//  - G, an n/2 x n generator matrix
// Output:
//  - Gi, an n/2 x n/2 invertible matrix composed of a random subset of the columns of G
void get_random_iset(mzd_t* G, mzd_t* Gis, rci_t* indices);

#endif
