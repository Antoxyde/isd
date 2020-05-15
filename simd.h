#ifndef ISD_POPCNT_H
#define ISD_POPCNT_H


#if defined(__AVX512F__) && defined(__AVX512BW__)

#include <m4ri.h>
#include <immintrin.h>


// When AX512 is enabled, we assume that the instance size is n = 1280
static inline uint16_t popcnt1280(void const* row) {

    uint16_t result = 0;

    __m512i m1 = _mm512_set1_epi8(0x55);
    __m512i m2 = _mm512_set1_epi8(0x33);
    __m512i m4 = _mm512_set1_epi8(0x0F);


    // first one
    __m512i a = _mm512_load_epi64(row);

    __m512i t1 = _mm512_sub_epi8(a, (_mm512_srli_epi16(a, 1) & m1));
    __m512i t2 = _mm512_add_epi8(t1 & m2, (_mm512_srli_epi16(t1, 2) & m2));
    __m512i t3 = _mm512_add_epi8(t2, _mm512_srli_epi16(t2, 4)) & m4;

    __m512i res_a =  _mm512_sad_epu8(t3, _mm512_setzero_si512());
    uint64_t* pca = (uint64_t*) &res_a;
    result += pca[0] + pca[1] + pca[2] + pca[3] +
         pca[4] + pca[5] + pca[6] + pca[7];

    // second one
    __m512i b = _mm512_load_epi64( (void*)(row + 64) );

    t1 = _mm512_sub_epi8(b, (_mm512_srli_epi16(b, 1) & m1));
    t2 = _mm512_add_epi8(t1 & m2, (_mm512_srli_epi16(t1, 2) & m2));
    t3 = _mm512_add_epi8(t2, _mm512_srli_epi16(t2, 4)) & m4;

    __m512i res_b =  _mm512_sad_epu8(t3, _mm512_setzero_si512());
    uint64_t* pcb = (uint64_t*) &res_b;
    result += pcb[0] + pcb[1] + pcb[2] + pcb[3] +
         pcb[4] + pcb[5] + pcb[6] + pcb[7];


    // third one
    __m256i c = _mm256_loadu_si256( (void*)(row + 128) );

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
    __m256i lo = _mm256_and_si256(c, low_mask);
    __m256i hi = _mm256_and_si256(_mm256_srli_epi16(c, 4), low_mask);
    __m256i popcnt1 = _mm256_shuffle_epi8(lookup1, lo);
    __m256i popcnt2 = _mm256_shuffle_epi8(lookup2, hi);
    __m256i res_c = _mm256_sad_epu8(popcnt1, popcnt2);
    uint64_t* pcc = (uint64_t*) &res_c;

    result += pcc[0] + pcc[1] + pcc[2] + pcc[3];

    return result;
}

static inline int cl_and_popcnt1280(mzd_t const* M1, mzd_t const* M2, rci_t i, rci_t j) {

    void* r1 = mzd_row(M1, i);
    void* r2 = mzd_row(M2, j);

    __m512i a1 = _mm512_load_epi64( r1 );
    __m512i a2 = _mm512_load_epi64( r2 );
    __m512i a = _mm512_xor_epi64(a1, a2);

    __m512i b1 = _mm512_load_epi64( (void*)(r1 + 64 ) );
    __m512i b2 = _mm512_load_epi64( (void*)(r2 + 64 ) );
    __m512i b = _mm512_xor_epi64(b1, b2);

    __m256i c1 = _mm256_load_si256( (void*)(r1 + 128) );
    __m256i c2 = _mm256_load_si256( (void*)(r2 + 128) );
    __m256i c = _mm256_xor_si256(c1, c2);

    // here is a slighty modified copy paste of the above popcnt1280, to make sure its inlined
    uint16_t result = 0;

    __m512i m1 = _mm512_set1_epi8(0x55);
    __m512i m2 = _mm512_set1_epi8(0x33);
    __m512i m4 = _mm512_set1_epi8(0x0F);


    // first one
    __m512i t1 = _mm512_sub_epi8(a, (_mm512_srli_epi16(a, 1) & m1));
    __m512i t2 = _mm512_add_epi8(t1 & m2, (_mm512_srli_epi16(t1, 2) & m2));
    __m512i t3 = _mm512_add_epi8(t2, _mm512_srli_epi16(t2, 4)) & m4;

    __m512i res_a =  _mm512_sad_epu8(t3, _mm512_setzero_si512());
    uint64_t* pca = (uint64_t*) &res_a;
    result += pca[0] + pca[1] + pca[2] + pca[3] +
         pca[4] + pca[5] + pca[6] + pca[7];

    // second one
    t1 = _mm512_sub_epi8(b, (_mm512_srli_epi16(b, 1) & m1));
    t2 = _mm512_add_epi8(t1 & m2, (_mm512_srli_epi16(t1, 2) & m2));
    t3 = _mm512_add_epi8(t2, _mm512_srli_epi16(t2, 4)) & m4;

    __m512i res_b =  _mm512_sad_epu8(t3, _mm512_setzero_si512());
    uint64_t* pcb = (uint64_t*) &res_b;
    result += pcb[0] + pcb[1] + pcb[2] + pcb[3] +
         pcb[4] + pcb[5] + pcb[6] + pcb[7];


    // third one
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
    __m256i lo = _mm256_and_si256(c, low_mask);
    __m256i hi = _mm256_and_si256(_mm256_srli_epi16(c, 4), low_mask);
    __m256i popcnt1 = _mm256_shuffle_epi8(lookup1, lo);
    __m256i popcnt2 = _mm256_shuffle_epi8(lookup2, hi);
    __m256i res_c = _mm256_sad_epu8(popcnt1, popcnt2);
    uint64_t* pcc = (uint64_t*) &res_c;

    result += pcc[0] + pcc[1] + pcc[2] + pcc[3];

    return result;
}

#else

static inline uint16_t popcnt(void* data, int size) {

    int i = 0;
    uint16_t result = 0;
    uint64_t* values_64 = (uint64_t*) data;
    int len = size / 8;
    int remaining = size % 8;

    for (i = 0; i < len; i++)
        result += __builtin_popcountll(values_64[i]);

    uint8_t* values_8 = (uint8_t*) (data + size * 8);

    for (i = 0; i < remaining; i++)
        result += __builtin_popcount(values_8[i]);

    return result;
}

#endif

#endif
