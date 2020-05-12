#ifndef ISD_UTILS_H
#define ISD_UTILS_H

#include <m4ri.h>

// Load the challenge given by the filename
// return 1 if everything went fine, 0 otherwise
int load_challenge(char* filename, mzd_t* G, mzd_t* H);

#endif



