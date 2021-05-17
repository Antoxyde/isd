#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <m4ri/m4ri.h>
#include "utils.h"

void test_precomp(mzd_t* H, uint64_t* H1tab, uint64_t* H2tab, uint64_t n, uint64_t k) {
    
    for (uint64_t val = 0; val < (1ULL << k); val++) {
        uint64_t result = 0;
        mzd_t* vec = mzd_init(1,n);
        mzd_t* res = mzd_init(1,n-k);

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

        result = H2tab[val >> n/2] ^ H1tab[val & ((1 << n/2) - 1)];

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
    
    uint64_t n = 40, k = 20;
    mzd_t* M = mzd_init(k, n-k);
    mzd_t* M2 = mzd_init(k, n-k);
    mzd_t* Ink = mzd_init(n-k,n-k);

    for (int i = 0; i < n-k; i++) mzd_write_bit(Ink, i, i, 1);
    
    do {
        mzd_randomize(M);
        mzd_copy(M2, M);
    } while (mzd_echelonize(M2, 0) != k);

    mzd_t* H = mzd_concat(NULL,  M, Ink); // H = [M | I_{n-k}]
    mzd_t* Ht = mzd_transpose(NULL, H);

    mzd_free(M2);

    uint64_t *Ht1tab = malloc(sizeof(uint64_t) * (1ULL << n/2));
    uint64_t *Ht2tab = malloc(sizeof(uint64_t) * (1ULL << n/2));
    CHECK_MALLOC(Ht1tab); CHECK_MALLOC(Ht2tab);
    
    mzd_t* vec_high, *vec_low, *res;
    vec_high = mzd_init(1, n);
    vec_low = mzd_init(1, n);
    res = mzd_init(1, n-k);

    uint64_t* alias_low = mzd_row(vec_low,0);
    uint64_t* alias_high = mzd_row(vec_high,0);
    
    mzd_info(H, 0);
    mzd_info(Ht, 0);

    printf("H tabulation precomputation...\n");
    for (*alias_low = 0; *alias_low < (1ULL << n/2); (*alias_low)++) {

        *alias_high = *alias_low << n/2;
        mzd_mul(res, vec_low, Ht, 0);
        Ht1tab[*alias_low] = *mzd_row(res, 0);

        mzd_mul(res, vec_high, Ht, 0);
        Ht2tab[*alias_low] = *mzd_row(res, 0);
    }
    printf("H tabulation precomputation done.\n");

    uint64_t t, old;

    //test_precomp(Ht, Ht1tab, Ht2tab, n, k);

    // MAP s -> lowest e s.t. s = eH
    uint64_t* SM = calloc(1ULL << (n-k), sizeof(uint64_t));
    CHECK_MALLOC(SM);

    uint64_t count = 0, total = 0;

    // Create all linear combinations of i rows of H, with i=1..n-k-1
    for (uint64_t i = 1; i < n; i++) { 

        uint64_t w = (1ULL << i) - 1; // i first bits to 1
        uint64_t syn = Ht2tab[w >> n/2] ^ Ht1tab[w & ((1 << n/2) - 1)];
        
        uint64_t nbiter = binomial(n, i);
        printf("i = %ld, nbiter = %lu\n", i, nbiter);
        for (uint64_t j = 0 ; j < nbiter; j++) {
            total++; 
            if (SM[syn] == 0) {
                count++;
                SM[syn] = w;
            }

            if (count == (1ULL << (n-k))) { 
                break;
            }
            
            if (j < nbiter-1) {
                old = w;
                // get next permutation of i bits, see https://graphics.stanford.edu/~seander/bithacks.html#NextBitPermutation
                t = w | (w - 1);
                w = (t + 1) | (((~t & -~t) - 1) >> (__builtin_ctz(old) + 1));  
                
                // Get the 2 indices that changed from last iteration
                uint64_t changes = old ^ w;
                // trailing zeroes, starting at LSB
                int i1 = __builtin_ctz(changes);
                // trailing zeroes, starting at MSB
                int i2 = 31 - __builtin_clz(changes);

                syn ^= (*mzd_row(Ht, i1)) ^ (*mzd_row(Ht, i2));
            }
        }

        if (count == (1ULL << (n-k))) { 
            break;
        }
    }

    if (count != (1ULL << (n-k))) {
        printf("Exhausted 2^n vectors and still have not found all syndromes O_o. Count is %lu\n", count);
    } else {
        printf("Syndrome table generated correctly! %lu total error vector for %lu syndromes\n", total, count);
    }




    
    mzd_free(H);
    mzd_free(Ht);
    free(Ht1tab);
    free(Ht2tab);
    free(SM);

    return 0;
}

