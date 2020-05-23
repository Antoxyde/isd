#ifndef ISD_ISET_H
#define ISD_ISET_H

#include <m4ri.h>

// Input:
//  - G, an n/2 x n generator matrix
//  - Gis, a n/2 x n/2 matrix, in which we will store the random iset
//  - indices, a set of values {0..n} to be shuffled
void get_random_iset(const mzd_t* Gt, mzd_t* Gis, mzd_t* Gist, rci_t* indices);

// Input:
//  - Glw, an n/2 x n matrix in systematic form
//  - indices, the permutation list to update
//  - affected_rows, a list to store all the rows that changed
void canteaut_next_iset(mzd_t* Glw, rci_t* perms, rci_t* affected_rows);

// Input:
//  - Glw, an n/2 x n/2 matrix, the redundant part of the current matrix
//  - indices, the permutation list to update
//  - affected_rows, a list to store all the rows that changed
void canteaut_next_iset_test(mzd_t* Glw, rci_t* perms, rci_t* affected_rows);

#endif
