#include <stdio.h>
#include <immintrin.h>
#include <stdint.h>

static inline uint64_t popcnt256(__m256i v)
{
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
  uint64_t* vals = (uint64_t*) res;

  return vals[0] + vals[1] + vals[2] + vals[3]
}


int main(void) {

    uint64_t a,b,c,d;
    a = 3;
    b = 3;
    c = 3;
    d = 3;

    __m256i mm = _mm256_loadu_si256(&a);
    __m256i r = popcnt256(mm);

    uint64_t* res = (uint64_t*)&r;
    int cnt = res[0] + res[1] + res[2] + res[3];
    printf("Res : %d\n", cnt);


    return 0;
}

