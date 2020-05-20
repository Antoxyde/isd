#include "utils.h"

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

    mzd_t* M = mzd_init(NROWS, NCOLS);
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

            mzd_write_bit(M, i, j, line[j] == '1' ? 1 : 0);
        }
    }

    mzd_concat(G, M, Ik); // G = [M| Ik]
    mzd_echelonize(G, 1);

    mzd_t* Mt = mzd_transpose(NULL, M);
    mzd_concat(H, Ik, Mt); // H = [Ik | M^t]

    fclose(fp);

    return 1;
}
