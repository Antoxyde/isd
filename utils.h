#ifndef ISD_UTILS_H
#define ISD_UTILS_H

#include <m4ri.h>

#define CHECK_MALLOC(X) if ((X) == NULL) {\
	fprintf(stderr, "Malloc failed at line %d in func %s\n", __LINE__, __func__);\
	exit(1);\
}



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



