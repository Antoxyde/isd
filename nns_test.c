#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <m4ri/m4ri.h>
//#include <malloc.h>

#include "utils.h"
#include "libpopcnt.h"
#include "xoshiro256starstar.h"
#include "buckets.h"

typedef struct TabulatedMat_ {
    uint64_t* t1,*t2;
} TabulatedMat;


/* w : vector
// Mtab : TabulatedMat*
// n : dimension of Mtab
// Compute wM
*/
#define TABULATED_MATMUL(w, Mtab, n) Mtab->t2[w >> n/2] ^ Mtab->t1[w & ((1 << n/2) - 1)]

int test_matrix_tabulation(mzd_t* Ht, TabulatedMat* Httab, uint64_t n, uint64_t k) {
    for (uint64_t val = 0; val < (1ULL << k); val++) {
        uint64_t result = 0;
        mzd_t* vec = mzd_init(1,n);
        mzd_t* res = mzd_init(1,n-k);

        uint64_t* alias_vec = mzd_row(vec, 0);
        *alias_vec = val;

        mzd_mul(res, vec, Ht, 0);
        uint64_t* alias_res = mzd_row(res, 0);

        result = TABULATED_MATMUL(val, Httab, n); 

        if (*alias_res != result) {
            printf("Matrix tabulation test failed for val = %lu\n", val);
            printf("alias_res : %lx\n", *alias_res);
            printf("result     : %lx\n", result);
            return 0;
        }
    }

    return 1;
}

int test_syndrome_table(uint64_t* SM,TabulatedMat* Httab, uint64_t n) {

    // Full test should be on all 2^n vectors, but depending on n it might be too expensive
    for (uint64_t e = 1; e < (1ULL << MIN(n, 26)); e++) { 
        uint64_t s = TABULATED_MATMUL(e, Httab, n); 
        
        uint64_t w1 = popcnt64_unrolled(&SM[s], 1);
        uint64_t w2 = popcnt64_unrolled(&e, 1);
        if (w1 > w2) {
            printf("Syndrome table test failed; for syndrome s=%lu SM[s] = %lu (wt : %lu), but a smaller error exists : %lu (wt : %lu) \n", s, SM[s], w1, e, w2);
            return 0;
        }
    }
    
    return 1;
}

/* Tabulate M , i.e. creates Mtab1 and Mtab2 s.t. xM = Mtab1[x LSBs] ^ Mtab2[x MSBs]*/
void tabulate_matrix(mzd_t* M, TabulatedMat* Mtab, uint64_t n, uint64_t k) {
    
    mzd_t* vec_high, *vec_low, *res;
    vec_high = mzd_init(1, n);
    vec_low = mzd_init(1, n);
    res = mzd_init(1, n-k);

    uint64_t* alias_low = mzd_row(vec_low,0);
    uint64_t* alias_high = mzd_row(vec_high,0);
    
    // For each possible n/2 LSBs and n/2 MSBs, precompute its matrix mult
    for (*alias_low = 0; *alias_low < (1ULL << n/2); (*alias_low)++) {

        *alias_high = *alias_low << n/2;
        mzd_mul(res, vec_low, M, 0);
        Mtab->t1[*alias_low] = *mzd_row(res, 0);

        mzd_mul(res, vec_high, M, 0);
        Mtab->t2[*alias_low] = *mzd_row(res, 0);
    }

}

/* Populates SM s.t. given a syndrome s, SM[s] = e where eH^t = s and e as the lowest weight possible */
int create_syndrome_table(uint64_t* SM, TabulatedMat* Httab, uint64_t n, uint64_t k) {
    
    uint64_t count = 0, total = 0, t, old, w, syn, nbiter, i, j;

    // Create all linear combinations of i rows of H, with i=1..n-k-1
    for (i = 1; i < n; i++) { 

        w = (1ULL << i) - 1; // i first bits to 1
        
        nbiter = binomial(n, i);
        for (j = 0 ; j < nbiter; j++) {
            total++; 
            
            // syn = wH^t
            syn = TABULATED_MATMUL(w, Httab, n); 
            if (SM[syn] == 0) {
                count++;
                SM[syn] = w;
            }

            if (count == (1ULL << (n-k))) { 
                break;
            }

            // get next permutation of i bits, see https://graphics.stanford.edu/~seander/bithacks.html#NextBitPermutation
            old = w; 
            t = w | (w - 1);
            w = (t + 1) | (((~t & -~t) - 1) >> (__builtin_ctzll(old) + 1));  
        }

        if (count == (1ULL << (n-k))) { 
            break;
        }
    }

    if (count != (1ULL << (n-k))) {
        printf("Exhausted 2^n vectors and still have not found all syndromes O_o. Count is %lu\n", count);
        return 0;
    } 
    
    return 1;
}

void make_stats(uint64_t* SM, TabulatedMat* Httab, mzd_t* G, uint64_t n, uint64_t k, uint64_t nbvec) {

    if (nbvec <= k) {
        printf("Nbvec has to be at least greater than k. (nbvec=%lu,k=%lu)\n", nbvec, k);
        return;
    }

    bucket** buckets = bucket_init(1ULL << k, 1ULL << (nbvec - k),32);
    uint64_t s, x, xc, e, w, i;
    uint64_t min_blen, max_blen;
    double av_blen;
    uint64_t perc = (1ULL << MIN(nbvec, n)) / 100;
    uint64_t *weights = calloc(sizeof(uint64_t) , k+1);
    
    printf("Populating buckets with 2^%lu vectors..\n", MIN(nbvec,n));

    for (i = 0; i < (1ULL << MIN(nbvec, n)); i++) {
        x = xoshiro256starstar_random() & ((1ULL << n) - 1);
        s = TABULATED_MATMUL(x, Httab, n);
        e = SM[s];
        xc = (x ^ e) & ((1ULL << k) - 1);
        bucket_put(buckets, xc, x);

        if (i % perc == 0) {
            printf("\r%lu%%", i / perc);
            fflush(stdout);
        }
    }
    
    printf("\r");
    fflush(stdout);

    mzd_t* cw = mzd_init(1, n);
    mzd_t* v = mzd_init(1, k);
    uint64_t* v_alias = mzd_row(v, 0);
    uint64_t* cw_alias = mzd_row(cw, 0);
    uint64_t dmin = n - k + 1; 
    uint64_t nb_empty_bucket = 0;
    double dinter = 0.0;
    perc = (1ULL << k) / 100;
    
    min_blen = max_blen = buckets[0]->curlen;
    av_blen = 0;
    
    printf("Making stats for each bucket..\n");
    // For each bucket
    for (i = 0; i < (1ULL << k); i++) {

        *v_alias = i;
        mzd_mul(cw, v, G, 0); // TODO : tabulate G also?
        w = popcnt64_unrolled(cw_alias, 1);
        if (w > 0) {
            dmin = MIN(dmin,  w);
        }

        bucket* b = buckets[i];
        double tot = 0.0;
        if (b && b->curlen > 1) {

            for (int j = 0; j < b->curlen; j++) {
                for (int k = j+1; k < b->curlen; k++) {
                    uint64_t r = b->tab[j] ^ b->tab[k];
                    uint64_t w = (double)popcnt64_unrolled(&r, 1);
                    tot += w;
                    weights[w] += 1;
                }
            }

            tot /= ((double)b->curlen * (b->curlen - 1))/2;
        } else {
            nb_empty_bucket++;
        }
        
        if (b) {
            av_blen += b->curlen;
            min_blen = MIN(min_blen, b->curlen);
            max_blen = MAX(min_blen, b->curlen);
        }

        dinter += tot;

        if (i % perc == 0) {
            printf("\r%lu%%", i / perc);
            fflush(stdout);
        }
    }

    printf("\r");
    fflush(stdout);
    
    av_blen /= (double)(1 << k);
    dinter /= (double)((1 << k) - nb_empty_bucket);
    
    printf("n                                  %lu\n", n);
    printf("k                                  %lu\n", k);
    printf("Number of vector in F^2_n hashed   2^%lu\n", MIN(nbvec, n));
    printf("Theoretical min weight (GV bound)  %lu\n", gv_bound(n, k));    
    printf("Practical min weight               %lu \n", dmin);
    printf("Average internal distance          %lf\n", dinter);
    printf("Number of empty bucket             %lu\n", nb_empty_bucket);
    printf("Min elements / bucket              %lu\n", nb_empty_bucket == 0 ? min_blen : 0);
    printf("Max elements / bucket              %lu\n", max_blen);
    printf("G                                  ");
    for (int i = 0; i < k; i++)
        printf("%lX ", *(uint64_t*)(mzd_row(G, i)));

    printf("\n");

    printf("Weight distrib                      ");
    for (int i = 0; i < k+1; i++) {
        printf("%d:%lu ", i, weights[i]);
    }
    printf("\n");

    
    bucket_free(buckets, 1ULL << k);
    mzd_free(cw);
    free(weights);
}

int main(int argc, char** argv) {
    
    //xoshiro256starstar_random_set((uint64_t[4]){1,3,3,8});
    if (argc < 4) {
        printf("Usage : %s n k nbvec\n" , argv[0]);
        return 1;
    }
    
    uint64_t n = atoll(argv[1]);
    uint64_t k = atoll(argv[2]);
    uint64_t nbvec = atoll(argv[3]);

    mzd_t* Ink = mzd_init(n-k,n-k);
    mzd_t* Ik = mzd_init(k,k);
    // n-k identity matrix
    for (int i = 0; i < n-k; i++) mzd_write_bit(Ink, i, i, 1);
    for (int i = 0; i < k; i++) mzd_write_bit(Ik, i, i, 1);
    
    mzd_t* M = get_random_fullrank(n-k, k);
    mzd_t* H = mzd_concat(NULL,  M, Ink); // H = [M | I_{n-k}]
    mzd_t* Ht = mzd_transpose(NULL, H); // H^t 
    mzd_t* Mt = mzd_transpose(NULL, M);
    mzd_t* G = mzd_concat(NULL, Ik, Mt);
    
    TabulatedMat Httab;
    Httab.t1 = malloc(sizeof(uint64_t) * (1ULL << n/2));
    Httab.t2 = malloc(sizeof(uint64_t) * (1ULL << n/2));
    CHECK_MALLOC(Httab.t1); CHECK_MALLOC(Httab.t2);

    // SM: s -> lowest e s.t. s = eH^t
    uint64_t* SM = calloc(1ULL << (n-k), sizeof(uint64_t));
    CHECK_MALLOC(SM);
    
    tabulate_matrix(Ht, &Httab, n, k);
    //if (!test_matrix_tabulation(Ht, &Httab, n, k)) return -1;

    create_syndrome_table(SM, &Httab, n, k);
    //if (!test_syndrome_table(SM, &Httab, n)) return -1;
    
    printf("Syndrome table generated.\n");

    make_stats(SM, &Httab, G, n, k, nbvec);

    mzd_free(H);
    mzd_free(Ht);
    mzd_free(G);
    mzd_free(Mt);
    free(Httab.t1);
    free(Httab.t2);
    free(SM);

    return 0;
}

