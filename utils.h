#ifndef ISD_UTILS_H
#define ISD_UTILS_H

#include <m4ri.h>

// Load the challenge given by the filename
// return H if loaded correctly, NULL otherwise
mzd_t* load_challenge(char* filename, mzd_t* H);

#endif



