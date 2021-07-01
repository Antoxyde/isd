#ifndef ISD_UTILS_H
#define ISD_UTILS_H

#include "m4ri/m4ri.h"

#define CHECK_MALLOC(X) do { \
        if ((X) == NULL) {\
	        fprintf(stderr, "Malloc failed at line %d in func %s\n", __LINE__, __func__);\
	        exit(1);\
        } \
    } while (0)


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


// Compute dst[i] ^= src[i] for i in 0..size
void mxor(uint64_t* dst, uint64_t* src, size_t size);

// Compute a^b but only on the hightest nbits, and cast the result in uint64
uint64_t uxor(uint64_t*a, uint64_t*b , size_t nbits);

void printbin(uint64_t* a, size_t nbits);

mzd_t* stern_reconstruct_cw(rci_t* min_comb, rci_t* column_perms, uint64_t* min_cw, uint64_t p);
mzd_t* prange_reconstruct_cw(rci_t row_min_cw, rci_t* column_perms, mzd_t* min_cw);

uint64_t binomial(uint64_t n, uint64_t k);

// Only for up to 64x64
void matrix_randomize(mzd_t* M, int r, int c);

// Only for up to 64x64
mzd_t* get_random_fullrank(int r, int c);
uint64_t gv_bound(uint64_t n, uint64_t k);

#endif



