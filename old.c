mzd_t* isd_lee_brickell(mzd_t* G, int niter) {

    int wt;
    rci_t n = G->ncols, i = 0;
    mzd_t* Gis = mzd_init(n/2, n/2);
    mzd_t* Gis_inv = mzd_init(n/2, n/2);
    mzd_t* Glw = mzd_init(n/2, n);
    mzd_t* Gt = mzd_transpose(NULL, G);
    mzd_t* Gist = mzd_init(n/2, n/2);

    int min_wt = n;
    mzd_t* min_cw = mzd_init(1, n);
    mzd_t* row_k = mzd_init(1, n);
    mzd_t* row_l = mzd_init(1, n);
    mzd_t* row_res = mzd_init(1,n);

    rci_t* indices = (rci_t*) malloc(sizeof(rci_t) * n);

    if (!indices) {
        fprintf(stderr, "Error in %s: failed to malloc %ld bytes.\n", __func__, sizeof(rci_t) * n);
        return NULL;
    }

    for (i = 0; i < n; i++) indices[i] = i;

    for (i = 0; i < niter; i++) {

        get_random_iset(Gt, Gis, Gist, indices); // Gis is n/2 x n/2
        mzd_inv_m4ri(Gis_inv, Gis, 0);
        mzd_mul(Glw, Gis_inv, G, 0); // Gi * G = Glw  ,  (n/2 x n/2) * (n/2 x n) => n/2 x n

        // Check all the CL of the rows for low cw
        for (rci_t k = 0; k < n/2; k++) {

            wt = popcnt(mzd_row(Glw, k), n/8 + (n % 8 != 0) );

            if (wt < min_wt) {
                printf("New min wt : %d\n", wt);
                min_wt = wt;
                mzd_copy_row(min_cw, 0, Glw, k);
            }


            for (rci_t l = 0; l < n/2; l++) {
                if (k != l) {

                    // TODO : optimized CL, without copying to a matrix
                    // eg. using 2 _mm512_loadu_epi8 and one

                    mzd_copy_row(row_l, 0, Glw, l);
                    mzd_copy_row(row_k, 0, Glw, k);
                    mzd_add(row_res, row_k, row_l);

                    wt = popcnt(mzd_first_row(row_res), n/8 + (n % 8 != 0) );

                    if (wt < min_wt) {
                        printf("New min wt : %d\n", wt);
                        min_wt = wt;
                        mzd_copy(min_cw, row_res);
                    }
                }
            }
        }

    }

    free(indices);

    mzd_free(row_l);
    mzd_free(row_k);
    mzd_free(row_res);

    mzd_free(Gt);
    mzd_free(Glw);
    mzd_free(Gist);
    mzd_free(Gis);
    mzd_free(Gis_inv);

    return min_cw;
}

mzd_t* isd_prange(mzd_t* G, int niter) {

    rci_t n = G->ncols, i = 0;
    mzd_t* Gis = mzd_init(n/2, n/2);
    mzd_t* Gis_inv = mzd_init(n/2, n/2);
    mzd_t* Glw = mzd_init(n/2, n);

    mzd_t* Gt = mzd_transpose(NULL, G);
    mzd_t* Gist = mzd_init(n/2, n/2);

    int min_wt = n;
    mzd_t* min_cw = mzd_init(1,n);

    rci_t* indices = (rci_t*) malloc(sizeof(rci_t) * n);

    if (!indices) {
        fprintf(stderr, "Error in %s: failed to malloc %ld bytes.\n", __func__, sizeof(rci_t) * n);
        return NULL;
    }

    for (i = 0; i < n; i++) indices[i] = i;

    for (i = 0; i < niter; i++) {

        get_random_iset(Gt, Gis, Gist, indices); // Gis is n/2 x n/2

        mzd_inv_m4ri(Gis_inv, Gis, 0);
        mzd_mul(Glw, Gis_inv, G, 0); // Gi * G = Glw  ,  (n/2 x n/2) * (n/2 x n) => n/2 x n

        // Check all the rows of Glw for low codewords
        for (rci_t j = 0; j < n/2; j++) {

            void* row = mzd_row(Glw, j);
            long int wt = popcnt(row, n/8 + (n % 8 != 0) );

            if (wt < min_wt) {
                printf("New min wt : %ld\n", wt);
                min_wt = wt;
                mzd_copy_row(min_cw, 0, Glw, j);
            }
        }

    }

    free(indices);

    mzd_free(Gist);
    mzd_free(Gt);
    mzd_free(Glw);
    mzd_free(Gis);
    mzd_free(Gis_inv);

    return min_cw;
}

mzd_t* isd_prange_canteaut_naive(mzd_t* G, int niter) {

    rci_t n = G->ncols, i = 0, j = 0;
    void* row = NULL;
    long int wt = 0;

    rci_t* column_perms_copy =  (rci_t*) malloc(sizeof(rci_t) * n);
    rci_t* column_perms = (rci_t*) malloc(sizeof(rci_t) * n);
    rci_t* affected_rows = (rci_t*) malloc(sizeof(rci_t) * n/2);

    if (!column_perms || !column_perms_copy || !affected_rows) {
        fprintf(stderr, "Error in %s:  malloc failed.\n", __func__);
        return NULL;
    }
    for (i = 0; i < n; i++) column_perms[i] = i;

    int min_wt = n;
    mzd_t* min_cw = mzd_init(1,n);
    mzd_t* Glw = mzd_copy(NULL, G);

    for (i = 0; i < niter; i++) {

        canteaut_next_iset(Glw, column_perms, affected_rows);

        // Check all the rows that changed since last iset for low codewords
        for (j = 0; affected_rows[j] >= 0; j++) {

            row = mzd_row(Glw, affected_rows[j]);
            //wt = 1 + popcnt(row + (n/16) /* n/2 bits, so n/16 bytes */ , n/16 + (n % 8 != 0) );
            wt = 1 + popcnt64_unrolled(row + (n/16), 10);

            if (wt < min_wt) {
                printf("New min wt : %ld\n", wt);
                min_wt = wt;

                // Save our new lowest row and all the permutations made until now
                mzd_copy_row(min_cw, 0, Glw, affected_rows[j]);
                memcpy(column_perms_copy, column_perms, n * sizeof(rci_t));
            }
        }
    }

    // Since we applied many columns permutations, which were all logged in perms
    // we have to permute back to get a valid codeword
    mzd_t* result = mzd_copy(NULL, min_cw);
    for (i = 0; i < n; i++) {
        if (i != column_perms_copy[i]) {
            printf("Deswap %d <-> %d\n", i, column_perms_copy[i]);
            mzd_write_bit(result, 0, column_perms_copy[i], mzd_read_bit(min_cw, 0, i));
        }
    }


    mzd_free(Glw);
    mzd_free(min_cw);

    free(affected_rows);
    free(column_perms);
    free(column_perms_copy);

    return result;
}


