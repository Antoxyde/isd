#ifndef ISD_UTILS_H
#define ISD_UTILS_H

#include <m4ri.h>

// Shuffle the size first entries of the values array
void fisher_yates_shuffle(rci_t* values, size_t size);

// Load the challenge given by the filename
// return 1 if everything went fine, 0 otherwise
int load_challenge(char* filename, mzd_t* G, mzd_t* H);

// Check whether the left squared part of M
// is the identity matrix or not
// ie if M is in systematic form
int left_is_identity(const mzd_t* M);

void print_cw(mzd_t* cw);

void rref_to_systematic(mzd_t* M, rci_t* perms);

#endif



