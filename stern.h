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

#define RADIX_WIDTH 2
#define RADIX_LEN (L/2)

#define STERN_GET(tab, index) ((tab[index >> 6]) >> (index & 0x3f)) & 1
#define STERN_SET_ONE(tab, index) (tab[index >> 6]) |= (1ULL << (index & 0x3f))

// Represent an LC array
typedef struct lc_tab_ {
    uint8_t p;
    size_t current_size;
    void* lcs;
} lc_tab;

// Represent a linear combination
typedef struct lc_ {
    uint32_t delta;
    uint16_t indexes[];
} lc;

#define LC_TAB_GET(tab, i) ((lc*)((tab)->lcs + (i) * (sizeof(uint32_t) + sizeof(uint16_t) * (tab)->p)))
#define LC_CALLOC(nelem, p) (lc*)(calloc(sizeof(uint32_t) + sizeof(uint16_t) * (p) , nelem))
#define LC_MALLOC(nelem, p) (lc*)(malloc((sizeof(uint32_t) + sizeof(uint16_t) * (p)) * nelem))

/* Input:
 *  - G, a n/2 x n generator matrix
 *  - 
 * Output:
 * - min_cw, the lowest codeword found
 */
 mzd_t* stern(mzd_t* G, uint64_t niter);
int compare_lc(const void* a, const void* b);
lc_tab* denomsort_r(lc_tab* T, lc_tab* Ts, int64_t Tlen, uint64_t width, uint64_t pos, uint32_t* aux);
lc_tab* radixsort(lc_tab* T, lc_tab* Ts, int64_t Tlen, uint32_t *aux);

#endif
