#include "libpopcnt.h"
#include "utils.h"
#include "xoshiro256starstar.h"

void fisher_yates_shuffle(rci_t* values, size_t size) {

    size_t i,j, tmp;

    for (i = size - 1; i > 0; i--) {
        j  = rand() % (i + 1);
        tmp = values[i];
        values[i] = values[j];
        values[j] = tmp;
    }
}

int load_challenge(char* filename, mzd_t* G, mzd_t* H) {

    FILE* fp = NULL;
    rci_t i = 0, j = 0;

    // size of the non-identity part of H
    const int NROWS = G->nrows;
    const int NCOLS = G->ncols / 2;

    char line[NROWS + 2]; // line are at most NROWS char + '\n' + '\0'

    fp = fopen(filename, "r");

    if (!fp) {
        fprintf(stderr, "Error in %s: failed to open the file %s for reading.\n", __func__, filename);
        return 0;
    }

    // Skip the 5 first lines (comments or constant stuff)
    char buf[1200];

    for (i = 0; i < 5; i++) {
        if (!fgets(buf, 1200, fp)) {
            fprintf(stderr, "Error in %s: failed to skip line (perhaps the file is incomplete?).\n", __func__);
            return 0;
        }
    }

    mzd_t* A = mzd_init(NROWS, NCOLS);
    mzd_t* Ik = mzd_init(NROWS, NCOLS);

    for (i = 0; i < NROWS; i++) {

        if (!fgets(line, NROWS + 2, fp)) {
            fprintf(stderr, "Error in %s: failed to read line (perhaps the file is incomplete?).\n", __func__);
            return 0;
        }

        for (j = 0; j < NCOLS; j++) {

            if (i == j) {
                mzd_write_bit(Ik, i, j, 1);
            }

            mzd_write_bit(A, i, j, line[j] == '1' ? 1 : 0);
        }
    }

    mzd_concat(G, A, Ik); // G = [A| Ik]
    mzd_echelonize(G, 1);

    mzd_t* Mt = mzd_transpose(NULL, A);
    mzd_concat(H, Ik, Mt); // H = [Ik | M^t]

    mzd_free(A);
    mzd_free(Ik);
    mzd_free(Mt);

    fclose(fp);

    return 1;
}

int left_is_identity(const mzd_t* A) {

    rci_t n = A->nrows;

    for (rci_t i = 0; i < n; i++) {
        if (mzd_read_bit(A, i, i) != 1 || popcnt(mzd_row(A, i), n/8) > 1) {
            return 0;
        }
    }

    return 1;
}

void print_cw(mzd_t* cw) {
    for (int i = 0; i < cw->ncols; i++) printf(mzd_read_bit(cw, 0, i) ? "1" : "0");
    printf("\n");
}

void rref_to_systematic(mzd_t* A, rci_t* perms) {

    rci_t n = A->ncols, cur, old;
    mzd_t* Mt = mzd_init(n, n/2);

    while (!left_is_identity(A)) {

        mzd_transpose(Mt, A);

        // Find a column with only one 1 in the right part of A
        for (cur = n/2; cur < n && popcnt(mzd_row(Mt, cur), n/16) != 1; cur++);

        // Find a column with more than one 1 in the left part of A
        for (old = 0; old < n/2 && popcnt(mzd_row(Mt, old), n/16) == 1; old++);

        // And swap them
        mzd_col_swap(A, old, cur);
        rci_t tmp = perms[cur];
        perms[cur] =  perms[old];
        perms[old] = tmp;

        mzd_echelonize(A, 1);
    }

    mzd_free(Mt);
}

void mxor(uint64_t* dst, uint64_t* src, size_t size) {
    for (size_t i = 0; i < size; i++)
        dst[i] ^= src[i];
}

uint64_t uxor(uint64_t*a, uint64_t*b , size_t nbits) {
    return (a[0] >> (64 - nbits)) ^ (b[0] >> (64 - nbits));
}

void printbin(uint64_t* a, size_t nbits) {

    for (size_t i = 0; i < nbits; i++) {
        printf("%ld", (a[i/64] >> ( 63 - (i%64))) & 1);
    }

    printf("\n");
}


mzd_t* stern_reconstruct_cw(rci_t* min_comb, rci_t* column_perms, uint64_t* min_cw, uint64_t psum) {

    mzd_t* ident = mzd_init(1, 640 /* k */);
    for (uint64_t i = 0; i < psum; i++)
        mzd_write_bit(ident, 0, min_comb[i], 1);

    mzd_t* min_cw_m = mzd_init(1, 640 /* k */);
    memcpy(mzd_first_row(min_cw_m), min_cw, 80);

    mzd_t* min_cw_full = mzd_concat(NULL, ident, min_cw_m);

    mzd_t* result = mzd_copy(NULL, min_cw_full);
    for (rci_t i = 0; i < 1280 /* n */; i++) {
        if (i != column_perms[i]) {
            mzd_write_bit(result, 0, column_perms[i], mzd_read_bit(min_cw_full, 0, i));
        }
    }

    mzd_free(min_cw_full);
    mzd_free(min_cw_m);
    mzd_free(ident);

    return result;
}

mzd_t* prange_reconstruct_cw(rci_t row_min_cw, rci_t* column_perms, mzd_t* min_cw) {

    // Since Glw contains only the redundant part of the matrix for the computation
    // we have to concat the identity part back again to get the correct codeword
    mzd_t* ident = mzd_init(1, 640 /* k */);

    // row_min_cw contain the row number from which min_cw has been taken
    mzd_write_bit(ident, 0, row_min_cw, 1);
    mzd_t* min_cw_full = mzd_concat(NULL, ident, min_cw);

    mzd_t* result = mzd_copy(NULL, min_cw_full);
    for (int i = 0; i < 1280 /* n */; i++) {
        if (i != column_perms[i]) {
            mzd_write_bit(result, 0, column_perms[i], mzd_read_bit(min_cw_full, 0, i));
        }
    }

    mzd_free(ident);
    mzd_free(min_cw_full);

    return result;
}

// Does not take care about overflows
uint64_t binomial(uint64_t n, uint64_t k) {
    if (k > n) return 0;
    if (k == 0 || k == n) return 1;

    k = MIN(k, n - k);
    uint64_t c = 1;
    
    for (uint64_t i = 0; i < k; i++) {
        c = c * (n - i) / (i + 1);
    } 
    return c;
}


// Must have c < 64
void matrix_randomize(mzd_t* A, int r, int c) {
    for(int i = 0; i < r; i++) {
        uint64_t rnd = xoshiro256starstar_random() & ((1ULL << c) - 1);
        uint64_t* alias_row = mzd_row(A, i);
        *alias_row = rnd;
    }
}

mzd_t* get_random_fullrank(int r, int c) {

    mzd_t* A = mzd_init(r, c);
    mzd_t* M2 = mzd_init(r, c);

    do {
        matrix_randomize(A, r, c);
        mzd_copy(M2, A);
    } while (mzd_echelonize(M2, 0) != MIN(r,c));
    
    mzd_free(M2);
    return A;
}


uint64_t gv_bound(uint64_t n, uint64_t k) {
    uint64_t max = 1ULL << (n-k);
    uint64_t d, s = 0;
    for (d = 1;s < max; d++) {
        s += binomial(n, d);
    }
    return d;
}

void canteaut_chabaud(mzd_t* Glw, rci_t* column_perms) {
    
    uint64_t tmp, j, lambda, mu;
    uint64_t* word = NULL;

    // Find lambda, mu s.t. Glw[lambda, mu] == 1
    uint64_t randdata =  xoshiro256starstar_random();
    lambda = randdata % 640 /* K */;
    randdata >>= 10 /* ceil(log_2(K)) for K=640 */;

    word = Glw->rows[lambda];

    mu = randdata % 10 /* K/64 */;
    randdata >>= 4 /* ceil(log_2(10)) */;

    // Assuming a whole row can't be zero
    while (word[mu] == 0) {
        mu = (mu + 1) % 10 /* K/64 */;
    }

    // We need to take a random 1 in this word.
    // Rotate the word by a random value and take the lowest one it has.
    j  = randdata % 64;
    uint64_t val = ((word[mu] << (64 - j)) | (word[mu] >> j));
    j = (j + _tzcnt_u64(val)) % 64;

#if defined(AVX512_ENABLED)
    uint64_t big_mu = mu, small_mu = j;
#endif
    mu = mu * 64 + j;

   // Log the column swapping
    tmp = column_perms[lambda];
    column_perms[lambda] = column_perms[mu + 640 /* K */];
    column_perms[mu + 640 /* K */] = tmp;

    // Clear the bit at the intersection of the lambda'th row and the mu'th column
    // so we don't have to rewrite a 1 in mu'th column everytime we add the lambda'th row
    mzd_write_bit(Glw, lambda, mu, 0);

#if defined(AVX512_ENABLED)
    void* row_lambda = Glw->rows[lambda];

    // A row contains K bits, so we use 512 + 128 bits registers
    __m512i rlambda1 = _mm512_loadu_si512(row_lambda);
    __m128i rlambda2 = _mm_loadu_si128(row_lambda + 64 /* 512/(8 * sizeof(void)) */);

    mask[big_mu] = (1ULL << small_mu);
    __m512i m1 = _mm512_loadu_si512(mask);
    __m128i m2 = _mm_loadu_si128(((void*)mask) + 64 /* 512/(8 * sizeof(void)) */);
#endif

    // Add the lambda'th row to every other row that have a 1 in the mu'th column
    for (j = 0; j < 640 /* K */; j++) {


#if defined(AVX512_ENABLED)
        void* row = Glw->rows[j];

        // Load the whole row
        __m512i rj1 = _mm512_loadu_si512(row);
        __m128i rj2 = _mm_loadu_si128(row + 64 /* 512/(8 * sizeof(void))  */);

        // Check whether there is a one in column mu using the mask
        if (j != lambda && (_mm512_test_epi64_mask(rj1, m1) >= 1 || _mm_test_epi64_mask(rj2, m2) >= 1)) {

            // Perform row addition
            _mm512_storeu_si512(row, _mm512_xor_si512(rlambda1, rj1));
            _mm_storeu_si128(row + 64 /* 512/(8 * sizeof(void)) */, _mm_xor_si128(rlambda2, rj2));
#else
        if (j != lambda && mzd_read_bit(Glw, j, mu) == 1) {
            mzd_row_add(Glw, lambda, j);
#endif
        }
    }

     // Unclear the bit we removed earlier
    mzd_write_bit(Glw, lambda, mu, 1);

#if defined(AVX512_ENABLED)
    // Clear the mask for the next iteration
    mask[big_mu] = 0;
#endif
}
