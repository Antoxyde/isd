#ifndef ISD_POPCNT_H
#define ISD_POPCNT_H


#ifdef __AVX512F__

#include <m4ri.h>
#include <immintrin.h>

static inline int popcnt512(__m512i val) {
    // TODO
}

static inline int popcnt256(__m256i val) {
    // TODO
}

// When AX512 is enabled, we assume that the instance size is n = 1280
static inline int popcnt1280(void const* row) {

    __m512i a = _mm512_load_epi64(row);
    __m512i b = _mm512_load_epi64( (void*)(row + 64) );
    __m256i c = _mm256_loadu_si256( (void*)(row + 128) );

    return 0;
}

static inline int cl_and_popcnt1280(mzd_t const* a, mzd_t const* b, rci_t i, rci_t j) {

    __m512i a1 = _mm512_load_epi64( mzd_row(a, i) );
    __m512i b1 = _mm512_load_epi64( mzd_row(b, j) );
    __m512i c1 = _mm512_xor_epi64(a1, b1);

    __m512i a2 = _mm512_load_epi64( (void*)(mzd_row(a, i) + 64 ) );
    __m512i b2 = _mm512_load_epi64( (void*)(mzd_row(b, j) + 64 ) );
    __m512i c2 = _mm512_xor_epi64(a2, b2);

    __m256i a3 = _mm256_load_si256( (void*)(mzd_row(a, i) + 128) );
    __m256i b3 = _mm256_load_si256( (void*)(mzd_row(b, i) + 128) );
    __m256i c3 = _mm256_xor_si256(a3, b3);

    return 0;
}

#else

static inline int popcnt(void* data, int size) {
    // TODO
    return 0;
}

#endif

#endif
