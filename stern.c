#include "stern.h"
#include "table.h"
#include "utils.h"

mzd_t* isd_stern_canteaut_chabaud_p2(mzd_t* G, uint64_t niter, uint64_t sigma) {

    uint64_t p = 2;
    rci_t n = G->ncols, comb1[2], comb2[2];


    uint64_t* linear_comb = (uint64_t*)malloc(sizeof(uint64_t) * 10);
    rci_t* column_perms_copy =  (rci_t*) malloc(sizeof(rci_t) * n);
    rci_t* column_perms = (rci_t*) malloc(sizeof(rci_t) * n);
    CHECK_MALLOC(colmun_perms);
    CHECK_MALLOC(column_perms_copy);
    CHECK_MALLOC(linear_comb);

    for (i = 0; i < n; i++) column_perms[i] = i;

    mzd_t* Gtemp = mzd_copy(NULL, G);

    // Ensure that we work with a systematic generator matrix
    rref_to_systematic(Gtemp, column_perms);
    mzd_t* Glw = mzd_submatrix(NULL, Gtemp, 0, n/2, n/2, n);

    // We store v=(a,b,c),k=(Z[a]+Z[b]+Z[c])[:sigma]
    table* tab = table_init(1 << sigma, 1);

    uint64_t delta1, delta2;
    int min_wt = n;

    for (comb1[0] = 0; comb1[0]  < n/4; comb1[0]++) {
        for (comb1[1] = 0; comb1[1] < n/4; comb1[1]++) {
            if (comb1[0] != comb1[1]) {
                delta1 = mxor(mzd_row(Glw, comb1[0]), mzd_row(Glw, comb1[1]), sigma);
                table_insert(tab, &comb1, p, delta1);
            }
        }
    }

    for (comb2[0] = n/4; comb2[0]  < n/2; comb2[0]++) {
        for (comb2[1] = n/4; comb2[1] < n/2; comb2[1]++) {
            if (comb2[0] != comb2[1]) {

                delta2 = mxor(mzd_row(Glw, comb2[0]), mzd_row(Glw, comb2[1]), sigma);
                bucket* buck = table_retrieve_bucket(ht, delta2);

                if (buck) {
                    for (i = 0; i < buck->len; i++) {
                        if (buck->elems[i]->key == delta2) { // If we have a collision

#if defined(AVX512_ENABLED)
                            // TODO cl av512
#else

                            comb1 = (rci_t*)(buck->elems[i]->data);
                            mxor(linear_comb,linear_comb, 80);
                            mxor(linear_comb,mzd_row(Glw, comb1[0]), 80);
                            mxor(linear_comb, mzd_row(Glw, comb1[1]), 80);
                            mxor(linear_comb, mzd_row(Glw, comb2[0]), 80);
                            mxor(linear_comb, mzd_row(Glw, comb2[1]), 80);
#endif
                            // TODO DBG : check linear_comb[:sigma] est bien a zéro
                            int wt = 2*p + popcnt(linear_comb, 10);

                            if (wt < min_wt) {
                                // Save the new min weight and the indexes of th e linear combination to obtain it
                                min_wt = wt;
                                memcpy(min_comb, comb1, p*sizeof(rci_t));
                                mempcy(min_comb + p, comb2, p*sizeof(rci_t));
                                // TODO save la permutation des colonnes actuels
                            }
                        }
                    }
                }
            }
        }
    }


    // TODO: next iset using canteaut&chabaud

    // TODO: permute back the columns
}
