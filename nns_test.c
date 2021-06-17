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

/* Verify that tabulated matrix-vector product is okay
 * for 2^k vectors */
int test_matrix_tabulation(mzd_t* Ht, TabulatedMat* Httab, uint64_t n, uint64_t k) {

#ifdef DEBUG
    printf("Started test matrix tabulation..\n");
#endif
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

#ifdef DEBUG
    printf("Test matrix tabulation ok.\n");
#endif
    return 1;
}

/* Verify that syndrome table is correct
 * by looking at the 2^MIN(n, 26) first vectors of F_2^n */
int test_syndrome_table(uint64_t* SM,TabulatedMat* Httab, uint64_t n) {

#ifdef DEBUG
    printf("Started test syndrome table..\n");
#endif
    // Full test should be on all 2^n vectors, but depending on n it might be too expensive
    for (uint64_t e = 0; e < (1ULL << MIN(n, 26)); e++) { 
        uint64_t s = TABULATED_MATMUL(e, Httab, n); 
        
        uint64_t w1 = popcnt64_unrolled(&SM[s], 1);
        uint64_t w2 = popcnt64_unrolled(&e, 1);

        if (w1 > w2) {
            printf("Syndrome table test failed; for syndrome s=%lu SM[s] = %lu (wt : %lu), but a smaller error exists : %lu (wt : %lu) \n", s, SM[s], w1, e, w2);
            return 0;
        }
    }
    
#ifdef DEBUG
    printf("Test syndrome table ok.\n");
#endif
    return 1;
}

/* Tabulate M , i.e. creates Mtab1 and Mtab2 s.t. xM = Mtab1[x LSBs] ^ Mtab2[x MSBs]*/
void tabulate_matrix(mzd_t* M, TabulatedMat* Mtab, uint64_t n, uint64_t k) {
    
#ifdef DEBUG
    printf("Started tabulate matrix..\n");
#endif

    mzd_t* vec_high, *vec_low, *res;
    vec_high = mzd_init(1, n);
    vec_low = mzd_init(1, n);
    res = mzd_init(1, n-k);

    uint64_t* alias_low = mzd_row(vec_low,0);
    uint64_t* alias_high = mzd_row(vec_high,0);
    
    // For each possible n/2 LSBs and n/2 MSBs, precompute its matrix mult
    for (*alias_low = 0; *alias_low < (1ULL << (n - n/2)); (*alias_low)++) {
        
        if (*alias_low < (1ULL << n/2)) {
            mzd_mul(res, vec_low, M, 0);
            Mtab->t1[*alias_low] = *mzd_row(res, 0);
        }

        *alias_high = *alias_low << n/2;
        mzd_mul(res, vec_high, M, 0);
        Mtab->t2[*alias_low] = *mzd_row(res, 0);
    }

#ifdef DEBUG
    printf("Tabulate matrix ok.\n");
#endif
}

/* Populates SM s.t. given a syndrome s, SM[s] = e where eH^t = s and e as the lowest weight possible */
int create_syndrome_table(uint64_t* SM, TabulatedMat* Httab, uint64_t n, uint64_t k) {
    
#ifdef DEBUG
    printf("Started syndrome table creation..\n");
#endif
    uint64_t count = 1, total = 0, t, old, w, syn, nbiter, i, j;

    // Create all linear combinations of i rows of H, with i=1..n-k-1
    for (i = 1; i < n; i++) { 

        w = (1ULL << i) - 1; // i first bits to 1
        
        nbiter = binomial(n, i);
        for (j = 0 ; j < nbiter; j++) {
            total++; 
            
            // syn = wH^t
            syn = TABULATED_MATMUL(w, Httab, n); 
            if (syn != 0 && SM[syn] == 0) {
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
    
#ifdef DEBUG
    printf("Syndrome table creation ok.\n");
#endif
    return 1;
}

/* Compute various stats of the code represented by G
 * - GV-bound and real min distance
 * - intra-bucket distance (average, using random vectors so with potential full-collisions) */
void make_stats(uint64_t* SM, TabulatedMat* Httab, mzd_t* G, uint64_t n, uint64_t k, uint64_t nbvec) {

    if (nbvec <= k) {
        printf("Nbvec has to be at least greater than k. (nbvec=%lu,k=%lu)\n", nbvec, k);
        return;
    }

    bucket** buckets = bucket_init(1ULL << k, 1ULL << (nbvec - k),32);
    uint64_t s, x, key, e, w, i;
    uint64_t min_blen, max_blen;
    double av_blen;
    uint64_t *weights = calloc(sizeof(uint64_t) , n+1);
    
#ifdef PROGRESS
    uint64_t perc = (1ULL << nbvec) / 100;
    printf("Populating buckets with 2^%lu vectors..\n", nbvec);
#endif

    for (i = 1; i < (1ULL << nbvec); i++) {
        x = xoshiro256starstar_random() & ((1ULL << n) - 1);
        //x = i & ((1ULL << n) - 1);
        s = TABULATED_MATMUL(x, Httab, n);
        e = SM[s];
        key = (x ^ e) & ((1ULL << k) - 1);
        bucket_put(buckets, key, x);

#ifdef PROGRESS
        if (perc && i % perc == 0) {
            printf("\r%lu%%", i / perc);
            fflush(stdout);
        }
#endif

    }

#ifdef PROGRESS
    printf("\r");
    fflush(stdout);
    perc = (1ULL << k) / 100;
    printf("Making stats for each bucket..\n");
#endif

    mzd_t* cw = mzd_init(1, n);
    mzd_t* v = mzd_init(1, k);
    uint64_t* v_alias = mzd_row(v, 0);
    uint64_t* cw_alias = mzd_row(cw, 0);
    uint64_t dmin = n - k + 1; 
    uint64_t nb_empty_bucket = 0;
    double dinter = 0.0;
    
    min_blen = max_blen = buckets[0]->curlen;
    av_blen = 0;

    // For each bucket
    for (i = 0; i < (1ULL << k); i++) {
        
        // Compute the weight of all vG to get the minimal distance
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
                    weights[w] += 1;
                }
            }

        } else {
            nb_empty_bucket++;
        }
        
        if (b) {
            av_blen += b->curlen;
            min_blen = MIN(min_blen, b->curlen);
            max_blen = MAX(min_blen, b->curlen);
        }

        dinter += tot;
#ifdef PROGRESS 
        if (perc && i % perc == 0) {
            printf("\r%lu%%", i / perc);
            fflush(stdout);
        }
#endif
    }
    
#ifdef PROGRESS
    printf("\r");
    fflush(stdout);
#endif
    
    av_blen /= (double)(1 << k);
    dinter = 0.0;
    uint64_t tot = 0;
    for (i = 0; i < k; i++) {
        tot += weights[i];
        dinter += (weights[i] * i);
    }

    dinter /= tot;
    printf("Number of combinations             %lu\n", tot);
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
    for (int i = 0; i < n+1; i++) {
        printf("%lu ", weights[i]);
    }
    printf("\n");

    
    bucket_free(buckets, 1ULL << k);
    mzd_free(cw);
    free(weights);
}

/* Same as make_stats, but use the code `nbconcat` times in direct sum*/
void make_stats_directsum(uint64_t* SM, TabulatedMat* Httab, uint64_t np, uint64_t kp, uint64_t nbvec, uint64_t nbconcat) {
    
    uint64_t k = nbconcat * kp;
    uint64_t n = nbconcat * np;

    if (nbvec <= k) {
        printf("Nbvec has to be at least greater than k. (nbvec=%lu,k=%lu)\n", nbvec, k);
        return;
    }

    bucket** buckets = bucket_init(1ULL << k, 1ULL << (nbvec - k),32);
    uint64_t s, x, key, e, w, i, j;
    uint64_t min_blen, max_blen;
    double av_blen;
    uint64_t *weights = calloc(sizeof(uint64_t) , n+1);
    
#ifdef PROGRESS
    uint64_t perc = (1ULL << MIN(nbvec, n)) / 100;
    printf("Populating buckets with 2^%lu vectors..\n", MIN(nbvec,n));
#endif

    for (i = 0; i < (1ULL << MIN(nbvec, n)); i++) {
        x = xoshiro256starstar_random() & ((1ULL << n) - 1);
        s = 0;
        for (j = 0; j < nbconcat; j++) {
            s = TABULATED_MATMUL((x >> (j*kp)) & ((1ULL << kp) - 1), Httab, np);
            e |= SM[s] << k*kp;
        }
        key = (x ^ e) & ((1ULL << k) - 1);
        bucket_put(buckets, key, x);

#ifdef PROGRESS
        if (perc && i % perc == 0) {
            printf("\r%lu%%", i / perc);
            fflush(stdout);
        }
#endif

    }

#ifdef PROGRESS
    printf("\r");
    fflush(stdout);
    perc = (1ULL << k) / 100;
    printf("Making stats for each bucket..\n");
#endif

    mzd_t* cw = mzd_init(1, n);
    mzd_t* v = mzd_init(1, k);
    uint64_t* v_alias = mzd_row(v, 0);
    uint64_t* cw_alias = mzd_row(cw, 0);
    uint64_t nb_empty_bucket = 0;
    double dinter = 0.0;
    
    min_blen = max_blen = buckets[0]->curlen;
    av_blen = 0;

    // For each bucket
    for (i = 0; i < (1ULL << k); i++) {

        bucket* b = buckets[i];
        double tot = 0.0;
        if (b && b->curlen > 1) {
            //printf("bucket %lu, len %lu\n", i, b->curlen);
            for (int j = 0; j < b->curlen; j++) {
                for (int k = j+1; k < b->curlen; k++) {
                    uint64_t r = b->tab[j] ^ b->tab[k];

                    //printf("r:%08lX, j:%08lX, k: %08lX\n", r, b->tab[j], b->tab[k]);
                    uint64_t w = (double)popcnt64_unrolled(&r, 1);
                    tot += w;
                    printf("w=%lu,k+1=%u\n", w, k+1);
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
#ifdef PROGRESS 
        if (perc && i % perc == 0) {
            printf("\r%lu%%", i / perc);
            fflush(stdout);
        }
#endif
    }


#ifdef PROGRESS
    printf("\r");
    fflush(stdout);
#endif
    
    av_blen /= (double)(1 << k);
    dinter /= (double)((1 << k) - nb_empty_bucket);
    
    printf("n                                  %lu\n", n);
    printf("k                                  %lu\n", k);
    printf("Number of vector in F^2_n hashed   2^%lu\n", nbvec);
    printf("Theoretical min weight (GV bound)  %lu\n", gv_bound(n, k));    
    printf("Average internal distance          %lf\n", dinter);
    printf("Number of empty bucket             %lu\n", nb_empty_bucket);
    printf("Min elements / bucket              %lu\n", nb_empty_bucket == 0 ? min_blen : 0);
    printf("Max elements / bucket              %lu\n", max_blen);

    printf("\n");

    printf("Weight distrib                      ");
    for (int i = 0; i < n+1; i++) {
        printf("%lu ", weights[i]);
    }
    printf("\n");

    
    bucket_free(buckets, 1ULL << k);
    mzd_free(cw);
    free(weights);
}

void get_hamming_3(mzd_t** G, mzd_t** H, uint64_t* n, uint64_t* k) {
    *k = 4;
    *n = 7;
    *G = mzd_from_str(*k, *n, "1000011010010100101100001111");
    *H = mzd_from_str(*n-*k, *n, "101010101100110001111");
}

void get_perfect_golay(mzd_t** G, mzd_t** H, uint64_t* n, uint64_t* k) {
    *n = 23;
    *k = 12;
    *G = mzd_from_str(*k, *n, "101011100011000000000000101011100011000000000000101011100011000000000000101011100011000000000000101011100011000000000000101011100011000000000000101011100011000000000000101011100011000000000000101011100011000000000000101011100011000000000000101011100011000000000000101011100011");
    *H = mzd_from_str( *n-*k,*n,"1000000000011111001001001000000000011111001001001000000001100011101100001000000001100011101100001000000110010001111000001000001001110101010000001000010110111100000000001000010110111100000000001000010110111100000000001000010110111100000000001111100100101");
}


int load_code(char** argv, int argc, mzd_t** G , mzd_t** H, uint64_t* n, uint64_t* k, uint64_t* nbvec) {

    if (argc != 3) {
        printf("Usage : %s load nbvec\n", argv[0]);
        return 0;
    } 

    *nbvec = atoll(argv[2]);
    get_hamming_3(G, H, n, k);
    return 1;
}

int random_code(char** argv, int argc, mzd_t** G , mzd_t** H, uint64_t* n, uint64_t* k, uint64_t* nbvec) {

    if (argc != 5) {

        printf("Usage : %s rand n k nbvec\n", argv[0]);
        return 0;
    }


    *n = atoll(argv[2]);
    *k = atoll(argv[3]);
    *nbvec = atoll(argv[4]);

    mzd_t* M = get_random_fullrank(*n-*k, *k);

    mzd_t* Ink = mzd_init(*n-*k,*n-*k);
    mzd_t* Ik = mzd_init(*k,*k);

    // n-k identity matrix
    for (int i = 0; i < *n-*k; i++) mzd_write_bit(Ink, i, i, 1);
    for (int i = 0; i < *k; i++) mzd_write_bit(Ik, i, i, 1);

    *H = mzd_concat(NULL,  M, Ink); // H = [M | I_{n-k}]
    mzd_t* Mt = mzd_transpose(NULL, M);
    *G = mzd_concat(NULL, Ik, Mt);

    mzd_free(Ik);
    mzd_free(Ink);
    mzd_free(M);
    mzd_free(Mt);

    return 1;
}

int main(int argc, char** argv) {
    
    // seeding not taken in account 
    //xoshiro256starstar_random_set((uint64_t[4]){1,3,3,8});
    
    mzd_t *G, *H;
    uint64_t n, k, nbvec;

    if (argc < 2) {
        printf("Usage : %s {rand, load}\n" , argv[0]);
        return 1;
    }

    if (!strcmp(argv[1], "rand")) {
        if (!random_code(argv, argc, &G, &H, &n, &k, &nbvec)) return -1;
    } else if (!strcmp(argv[1], "load")) {
        if (!load_code(argv, argc, &G, &H, &n, &k, &nbvec)) return -1;
    } else {
        printf("arg 1 as to be either 'rand' or 'load'.\n");
        return 2;
    }
    
    // From there, G, H, n, k, nbvec are initialised
    
    mzd_t* Ht = mzd_transpose(NULL, H);

    TabulatedMat Httab;
    Httab.t1 = malloc(sizeof(uint64_t) * (1ULL << n/2));
    Httab.t2 = malloc(sizeof(uint64_t) * (1ULL << (n - n/2)));
    CHECK_MALLOC(Httab.t1); CHECK_MALLOC(Httab.t2);

    // SM: s -> lowest e s.t. s = eH^t
    uint64_t* SM = calloc(1ULL << (n-k), sizeof(uint64_t));
    CHECK_MALLOC(SM);
    
    tabulate_matrix(Ht, &Httab, n, k);
    if (!test_matrix_tabulation(Ht, &Httab, n, k)) return -1;

    create_syndrome_table(SM, &Httab, n, k);
    if (!test_syndrome_table(SM, &Httab, n)) return -1;
    
    printf("Syndrome table generated.\n");
    
    make_stats(SM, &Httab, G, n, k, nbvec);
    //make_stats_directsum(SM, &Httab, n, k, nbvec, 2);


    mzd_free(H);
    mzd_free(Ht);
    mzd_free(G);
    free(Httab.t1);
    free(Httab.t2);
    free(SM);

    return 0;
}

