#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <m4ri/m4ri.h>
#include "utils.h"

void test_precomp(mzd_t* H, uint64_t* H1tab, uint64_t* H2tab, uint64_t n, uint64_t k) {
    
    for (uint64_t val = 0; val < (1ULL << n-k); val++) {
        uint64_t result = 0, val = 0x5;
        mzd_t* vec = mzd_init(1,n-k);
        mzd_t* res = mzd_init(1,n);

        uint64_t* alias_vec = mzd_row(vec, 0);
        *alias_vec = val;

        mzd_mul(res, vec, H, 0);
        uint64_t* alias_res = mzd_row(res, 0);
        /* 
        printf("H:\n");
        mzd_info(H, 0);
        mzd_print(H);

        printf("vec:\n");
        mzd_info(vec, 0);
        mzd_print(vec);

        printf("res : \n");
        mzd_info(res, 0);
        mzd_print(res);
        */

        result = H2tab[val >> (n-k)/2] ^ H1tab[val & ((1 << (n-k)/2) - 1)];

        if (*alias_res != result) {
            printf("test_precomp failed for val = %lu\n", val);
            printf("alias_res : %lx\n", *alias_res);
            printf("result     : %lx\n", result);
            return;
        }
    }

    printf("test_precomp ok\n");
}


int main(void) {
    
    uint64_t n = 10, k = 5;
    mzd_t* M = mzd_init(n-k, k);
    mzd_t* M2 = mzd_init(n-k, k);
    mzd_t* Ink = mzd_init(n-k,n-k);

    for (int i = 0; i < n-k; i++) mzd_write_bit(Ink, i, i, 1);
    
    do {
        mzd_randomize(M);
        mzd_copy(M2, M);
    } while (mzd_echelonize(M2, 0) != k);

    mzd_t* H = mzd_concat(NULL,  M, Ink); // H = [M | I_{n-k}]
    mzd_free(M2);

    uint64_t nk_dec = (n-k)/2;

    uint64_t *H1tab = malloc(sizeof(uint64_t) * (1ULL << nk_dec));
    uint64_t *H2tab = malloc(sizeof(uint64_t) * (1ULL << nk_dec));
    CHECK_MALLOC(H1tab); CHECK_MALLOC(H2tab);
    
    mzd_t* vec_high, *vec_low, *res;
    vec_high = mzd_init(1, n-k);
    vec_low = mzd_init(1, n-k);
    res = mzd_init(1, n);

    uint64_t* alias_low = mzd_row(vec_low,0);
    uint64_t* alias_high = mzd_row(vec_high,0);

    mzd_info(H, 0);

    printf("Hitab precomputation...\n");
    for (*alias_low = 0; *alias_low < (1ULL << nk_dec); (*alias_low)++) {

        *alias_high = *alias_low << nk_dec;
        mzd_mul(res, vec_low, H, 0);
        H1tab[*alias_low] = *mzd_row(res, 0);

        mzd_mul(res, vec_high, H, 0);
        H2tab[*alias_low] = *mzd_row(res, 0);
    }
    printf("Hitab precomputation done.\n");

    uint64_t t, old;

    test_precomp(H, H1tab, H2tab, n, k);

    /*
    // Create all linear combinations of i rows of H, with i=1..n-k-1
    for (uint64_t i = 1; i < n-k; i++) { 

        uint64_t w = (1ULL << i) - 1; // i first bits to 1
        uint64_t s = H1tab[w & ((1U << nk_dec) - 1)] ^ H2tab[w >> nk_dec];

        for (uint64_t j = 0 ; j < binomial(n-k,i); j++) {

            old = w;
            // get next permutation of i bits, see https://graphics.stanford.edu/~seander/bithacks.html#NextBitPermutation
            t = w | (w - 1);
            w = (t + 1) | (((~t & -~t) - 1) >> (__builtin_ctz(v) + 1));  
            
            // Get the 2 indices that changed from last iteration
            uint64_t changes = old ^ w;
            // trailing zeroes, starting at LSB
            int i1 = __builtin_ctz(changes);
            // trailing zeroes, starting at MSB
            int i2 = 31 - __builtin_clz(changes);

            // TODO syn ^= H.row(i) ^ H.row(i2)
            syn ^= ((*mzd_row(H, i)) >> (64 - nk_dec)) ^  ((*mzd_row(H, j)) >> (64 - nk_dec))
        }
    }
    */
    return 0;
}

