#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <m4ri/m4ri.h>

#include "utils.h"
#include "libpopcnt.h"

int test_matrix_tabulation(mzd_t* H, uint64_t* H1tab, uint64_t* H2tab, uint64_t n, uint64_t k) {
    printf("Testing matrix tabulation..\n");    
    for (uint64_t val = 0; val < (1ULL << k); val++) {
        uint64_t result = 0;
        mzd_t* vec = mzd_init(1,n);
        mzd_t* res = mzd_init(1,n-k);

        uint64_t* alias_vec = mzd_row(vec, 0);
        *alias_vec = val;

        mzd_mul(res, vec, H, 0);
        uint64_t* alias_res = mzd_row(res, 0);

        result = H2tab[val >> n/2] ^ H1tab[val & ((1 << n/2) - 1)];

        if (*alias_res != result) {
            printf("Matrix tabulation test failed for val = %lu\n", val);
            printf("alias_res : %lx\n", *alias_res);
            printf("result     : %lx\n", result);
            return 0;
        }
    }
    
    printf("Ok.\n");

    return 1;
}

int test_syndrome_table(uint64_t* SM,uint64_t* Httab1, uint64_t* Httab2, uint64_t n, uint64_t k) {
    printf("Starting syndrome table test...\n");

    // Full test should be on all 2^n vectors, but depending on n it might be too expensive
    for (uint64_t e = 1; e < (1ULL << MIN(n, 26)); e++) { 
        uint64_t s = Httab2[e >> n/2] ^ Httab1[e & ((1 << n/2) - 1)];

        if (popcnt64_unrolled(&SM[s], n) > popcnt64_unrolled(&e, n)) {
            printf("Syndrome table test failed; for syndrome s=%lx, a smaller error exists : %lx (SM[s] = %lx)\n", s, e, SM[s]);
            return 0;
        }
    }
    
    printf("Ok.\n");
    return 1;
}

/* Tabulate M s.t. xM = Mtab1[x LSBs] ^ Mtab2[x MSBs]*/
void tabulate_matrix(mzd_t* M, uint64_t* Mtab1, uint64_t* Mtab2, uint64_t n, uint64_t k) {
    
    printf("Tabulating matrix...\n");
    mzd_t* vec_high, *vec_low, *res;
    vec_high = mzd_init(1, n);
    vec_low = mzd_init(1, n);
    res = mzd_init(1, n-k);

    uint64_t* alias_low = mzd_row(vec_low,0);
    uint64_t* alias_high = mzd_row(vec_high,0);

    for (*alias_low = 0; *alias_low < (1ULL << n/2); (*alias_low)++) {

        *alias_high = *alias_low << n/2;
        mzd_mul(res, vec_low, M, 0);
        Mtab1[*alias_low] = *mzd_row(res, 0);

        mzd_mul(res, vec_high, M, 0);
        Mtab2[*alias_low] = *mzd_row(res, 0);
    }

    printf("Done.\n");
}

int create_syndrome_table(uint64_t* SM, mzd_t* Ht, uint64_t* Httab1, uint64_t* Httab2, uint64_t n, uint64_t k) {
    
    printf("Creating syndrome table...\n");
    uint64_t count = 0, total = 0, t, old, w, syn, nbiter, i, j, changes;

    // Create all linear combinations of i rows of H, with i=1..n-k-1
    for (i = 1; i < n; i++) { 

        w = (1ULL << i) - 1; // i first bits to 1
        syn = Httab2[w >> n/2] ^ Httab1[w & ((1 << n/2) - 1)];
        
        nbiter = binomial(n, i);
        printf("wt(e) = %ld, nb possibilities = %lu\n", i, nbiter);
        for (j = 0 ; j < nbiter; j++) {
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
                changes = old ^ w;
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
        return 0;
    } 
    
    printf("Done.\n");
    return 1;
}

mzd_t* get_random_fullrank(int r, int c) {

    mzd_t* M = mzd_init(r, c);
    mzd_t* M2 = mzd_init(r, c);

    do {
        mzd_randomize(M);
        mzd_copy(M2, M);
    } while (mzd_echelonize(M2, 0) != MIN(r,c));
    
    mzd_free(M2);
    return M;
}
int main(void) {
    
    uint64_t n = 40, k = 20;

    mzd_t* Ink = mzd_init(n-k,n-k);
    // n-k identity matrix
    for (int i = 0; i < n-k; i++) mzd_write_bit(Ink, i, i, 1);
    
    mzd_t* M = get_random_fullrank(k, n-k);
    mzd_t* H = mzd_concat(NULL,  M, Ink); // H = [M | I_{n-k}]
    mzd_t* Ht = mzd_transpose(NULL, H);

    uint64_t *Httab1 = malloc(sizeof(uint64_t) * (1ULL << n/2));
    uint64_t *Httab2 = malloc(sizeof(uint64_t) * (1ULL << n/2));
    CHECK_MALLOC(Httab1); CHECK_MALLOC(Httab2);

    // MAP s -> lowest e s.t. s = eH^t
    uint64_t* SM = calloc(1ULL << (n-k), sizeof(uint64_t));
    CHECK_MALLOC(SM);
    
    tabulate_matrix(Ht, Httab1, Httab2, n, k);
    test_matrix_tabulation(Ht, Httab1, Httab2, n, k);

    create_syndrome_table(SM, Ht, Httab1, Httab2, n, k);
    test_syndrome_table(SM, Httab1, Httab2, n, k);

    printf("Syndrome table generated, all okay.\n");

    mzd_free(H);
    mzd_free(Ht);
    free(Httab1);
    free(Httab2);
    free(SM);

    return 0;
}

