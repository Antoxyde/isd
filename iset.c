#include "iset.h"
#include "utils.h"

void get_random_iset(const mzd_t* Gt, mzd_t* Gis, mzd_t* Gist, rci_t* indices) {

    rci_t n = Gt->nrows;

    // M4ri have no function to copy a column from a matrix to another
    // so we transpose and copy rows instead

    // while our matrix is not invertible
    do {

        // we choose n/2 random columns
        fisher_yates_shuffle(indices, n);

        for (rci_t col = 0; col < n/2; col++) {
            mzd_copy_row(Gist, col, Gt, indices[col]);
        }

    } while (mzd_echelonize(Gist, 0) <= ((n/2) - (n/128)) );

    // Since we copied rows from Gt, we have to tranpose back again
    mzd_transpose(Gis, Gist);
}

void canteaut_next_iset(mzd_t* Glw, rci_t* perms, rci_t* affected_rows) {

    // Canteaut improvement to derive a new iset from a previous one

    int pos = 0; // next free position to write in affected_rows
    rci_t n = Glw->ncols, lambda, mu, i, tmp;

    do {
        lambda = rand() % (n/2);
        mu = (rand() % (n/2) ) + (n/2);
    } while (mzd_read_bit(Glw, lambda, mu) == 0);

    mzd_col_swap(Glw, lambda, mu);

    tmp = perms[lambda];
    perms[lambda] = perms[mu];
    perms[mu] = tmp;

    for (i = 0; i < n/2; i++) {
        if (i != lambda && mzd_read_bit(Glw, i, lambda) == 1) {
            affected_rows[pos++] = i;
            mzd_row_add(Glw, lambda, i);
        }
    }

    affected_rows[pos] = -1;
}

void canteaut_next_iset_test(mzd_t* Glw, rci_t* perms, rci_t* affected_rows) {

    // Canteaut improvement to derive a new iset from a previous one

    int current = 0;
    rci_t n = Glw->ncols * 2, lambda, mu, i, tmp;

    do {
        lambda = rand() % (n/2);
        mu = rand() % (n/2);
    } while (mzd_read_bit(Glw, lambda, mu) == 0);

    tmp = perms[lambda];
    perms[lambda] = perms[mu + (n/2)];
    perms[mu + (n/2)] = tmp;

    for (i = 0; i < n/2; i++) {
        if (i != lambda && mzd_read_bit(Glw, i, mu) == 1) {
            affected_rows[current++] = i;
            mzd_row_add(Glw, lambda, i);
            mzd_write_bit(Glw, i, mu, 1);
        }
    }

    affected_rows[current] = -1;
}

