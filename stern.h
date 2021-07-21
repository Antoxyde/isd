#ifndef ISD_STERN_H
#define ISD_STERN_H

#if defined(__AVX512DQ__) && defined(__AVX512F__) && defined(__AVX512VL__)
#define AVX512_ENABLED
#endif

#include "m4ri/m4ri.h"
#include <stdint.h>

#define M       1         // Number of window
#define P1      4         // Stern's P parameter
#define P2      3         // Stern's P parameter
#define N       1280      // Code length
#define K       640       // Code dimension
#define L       20        // Stern's L paramter, window size 
#define CC      320       // Number of canteaut-chabaud iterations

#if defined(L2_THRESHOLD)
#define L2_SIZE_THRESHOLD       (1ULL << 16)
#endif /* L2_THRESHOLD*/

#define RADIX_WIDTH 2
#define RADIX_LEN (L/2)

#define STERN_GET(tab, index) ((tab[index >> 6]) >> (index & 0x3f)) & 1
#define STERN_SET_ONE(tab, index) (tab[index >> 6]) |= (1ULL << (index & 0x3f))

/* We can't have an array of struct with dynamically-sized in C without using malloc (Actually we can have `flexible array members` but not arrays of them)
* Calling malloc for each of `PX choose K/2` instances is not a possibility
* So we use a small hack by keeping the number of elements lc->indexes in a parent struct (lc_tab->p)
* And convenients macros to use this weird structure
* Admit its not really 
*/

// Represent an LC array
typedef struct lc_tab_ {
    uint8_t p; // used to have a dynamic size for the elements of lcs
    size_t current_size; // Current number of elements in lcs
    void* lcs; // pointer to an array of lc*
} lc_tab;


// Represent a linear combination
typedef struct lc_ {
    uint32_t delta; // Value on the window for this combination
    uint16_t indexes[]; // Row indexes to add to obtain the combination. Normalised to [0,K/2].
} lc;

// Convenient macros for the use of the above structures
#define LC_TAB_GET(tab, i) ((lc*)((tab)->lcs + (i) * (sizeof(uint32_t) + sizeof(uint16_t) * (tab)->p)))
#define LC_CALLOC(nelem, p) (lc*)(calloc(sizeof(uint32_t) + sizeof(uint16_t) * (p) , nelem))
#define LC_MALLOC(nelem, p) (lc*)(malloc((sizeof(uint32_t) + sizeof(uint16_t) * (p)) * nelem))

/* Input:
 *  - G, a n/2 x n generator matrix
 *  - time_sec, number of seconds to run the algorithm
 * Output:
 * - min_cw, the lowest codeword found
 */
 mzd_t* stern(mzd_t* G, uint64_t time_sec);
lc_tab* denomsort_r(lc_tab* T, lc_tab* Ts, uint64_t pos, uint32_t* aux);
lc_tab* radixsort(lc_tab* T, lc_tab* Ts, uint32_t *aux);

#endif
