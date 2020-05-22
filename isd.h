#ifndef ISD_ISD_H
#define ISD_ISD_H

#include <m4ri.h>

// Input:
//  - G, an n/2 x n generator matrix
//  - Gis, a n/2 x n/2 matrix, in which we will store the random iset
//  - indices, a set of values {0..n} to be shuffled
void get_random_iset(const mzd_t* Gt, mzd_t* Gis, mzd_t* Gist, rci_t* indices);


// Input:
//  - Glw, an n/2 x n matrix in systematic form
//  - indices, the permutation list to update
void canteaut_next_iset_naive(mzd_t* Glw, rci_t* perms);

// Input:
//  - Glw, an n/2 x n matrix in systematic form
//  - indices, the permutation list to update
//  - affected_rows, a list to store all the rows that changed
void canteaut_next_iset(mzd_t* Glw, rci_t* perms, rci_t* affected_rows);

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

#endif
