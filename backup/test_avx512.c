#include <immintrin.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>

int main(void) {

        uint64_t a[8];
        uint64_t b[8];
        uint64_t c[8];

        for (int i = 0; i < 8; i++) {
                a[i] = i;
                b[i] = i;
                c[i] = i;
        }

        /*
         // Segfault
        __m512i* aa = (__m512i*)a;
        __m512i* bb = (__m512i*)b;
        *(__m512i*)c = _mm512_xor_si512(*aa, *bb);
        */

        /*
         // Fonctionne mais apparement lent
        __m512i aa = _mm512_loadu_si512(a);
        __m512i bb = _mm512_loadu_si512(b);
        *(__m512i*)c = _mm512_xor_si512(aa,bb);
        */

        // fonctionne et rapide
        __m512i aa = _mm512_loadu_si512(a);
        __m512i bb = _mm512_loadu_si512(b);
        _mm512_storeu_si512(&b, _mm512_xor_si512(aa,bb));

        for (int i = 0; i < 8; i++) printf("%ld,", b[i]);
        printf("\n");

        return 0;
}

