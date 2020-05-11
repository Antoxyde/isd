#include "utils.h"


mzd_t* load_challenge(char* filename, mzd_t* H) {

    FILE* fp = NULL;
    int i = 0, j = 0;
    const int NCOLS = 640;
    const int NROWS = 640;
    char line[NROWS + 2]; // line are at most NROWS char + '\n' + '\0'

    if (H->nrows != NROWS || H->ncols != 2*NCOLS) {
        fprintf(stderr, "Error in %s: H should be a %dx%d matrix, but its a %dx%d one.\n", __func__, NROWS, 2*NCOLS, H->nrows, H->ncols);
        return NULL;
    }

    fp = fopen(filename, "r");

    if (!fp) {
        fprintf(stderr, "Error in %s: failed to open the file %s for reading.\n", __func__, filename);
        return NULL;
    }


    // Skip the 5 first lines
    for (i = 0; i < 5; i++) {
        if (!fgets(line, NCOLS + 2, fp)) {
            fprintf(stderr, "Error in %s: failed to skip line (perhaps the file is incomplete?).\n", __func__);
            return NULL;
        }
    }

    // 640 next lines are the non-identity part of H
    for (i = 0; i < NCOLS; i++) {

        if (!fgets(line, NROWS + 2, fp)) {
            fprintf(stderr, "Error in %s: failed to read line (perhaps the file is incomplete?).\n", __func__);
            return NULL;
        }

        for (j = 0; j < NROWS; j++) {

            // The 640 fist columns are the identity
            if (i == j)
                mzd_write_bit(H, i, j, 1);

            //
            mzd_write_bit(H, j, NCOLS + i, line[j] == '1' ? 1 : 0);
        }
    }

    fclose(fp);

    return H;
}
