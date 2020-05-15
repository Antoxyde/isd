#ifndef ISD_POPCNT_H
#define ISD_POPCNT_H


#ifdef __AVX512F__

#include <m4ri.h>
#include <immintrin.h>

// Mula's algorithm
// adapted from libpopcnt
static inline uint16_t popcnt256(__m256i v) {
  __m256i lookup1 = _mm256_setr_epi8(
      4, 5, 5, 6, 5, 6, 6, 7,
      5, 6, 6, 7, 6, 7, 7, 8,
      4, 5, 5, 6, 5, 6, 6, 7,
      5, 6, 6, 7, 6, 7, 7, 8
  );

  __m256i lookup2 = _mm256_setr_epi8(
      4, 3, 3, 2, 3, 2, 2, 1,
      3, 2, 2, 1, 2, 1, 1, 0,
      4, 3, 3, 2, 3, 2, 2, 1,
      3, 2, 2, 1, 2, 1, 1, 0
  );

  __m256i low_mask = _mm256_set1_epi8(0x0f);
  __m256i lo = _mm256_and_si256(v, low_mask);
  __m256i hi = _mm256_and_si256(_mm256_srli_epi16(v, 4), low_mask);
  __m256i popcnt1 = _mm256_shuffle_epi8(lookup1, lo);
  __m256i popcnt2 = _mm256_shuffle_epi8(lookup2, hi);

  __m256i res = _mm256_sad_epu8(popcnt1, popcnt2);
  uint64_t* vals = (uint64_t*) &res;

  return vals[0] + vals[1] + vals[2] + vals[3];
}

static inline uint16_t popcnt512(__m512i v) {
  __m512i m1 = _mm512_set1_epi8(0x55);
  __m512i m2 = _mm512_set1_epi8(0x33);
  __m512i m4 = _mm512_set1_epi8(0x0F);
  __m512i t1 = _mm512_sub_epi8(v, (_mm512_srli_epi16(v, 1) & m1));
  __m512i t2 = _mm512_add_epi8(t1 & m2, (_mm512_srli_epi16(t1, 2) & m2));
  __m512i t3 = _mm512_add_epi8(t2, _mm512_srli_epi16(t2, 4)) & m4;

  __m512i res =  _mm512_sad_epu8(t3, _mm512_setzero_si512());
  uint64_t* vals = (uint64_t*) &res;

  return vals[0] + vals[1] + vals[2] + vals[3] +
         vals[4] + vals[5] + vals[6] + vals[7];
}


// When AX512 is enabled, we assume that the instance size is n = 1280
static inline uint16_t popcnt1280(void const* row) {

    __m512i a = _mm512_load_epi64(row);
    __m512i b = _mm512_load_epi64( (void*)(row + 64) );
    __m256i c = _mm256_loadu_si256( (void*)(row + 128) );

    return popcnt512(a) + popcnt512(b) + popcnt256(c);
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

    return popcnt512(c1) + popcnt512(c2) + popcnt256(c3);
}

#else

static inline int popcnt(void* data, int size) {
    // TODO
    return 0;
}

#endif

#endif
