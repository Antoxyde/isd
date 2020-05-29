#include "stern.h"
#include "hashtable.h"
#include "utils.h"

mzd_t* isd_stern_canteaut_chabaud(mzd_t* G, uint64_t niter) {

    // Assuming p=2 and sigma=18
    uint64_t sigma = 18;
    uint64_t p = 2;
    rci_t n = G->ncols, comb1[p], comb2[p];

    rci_t* column_perms_copy =  (rci_t*) malloc(sizeof(rci_t) * n);
    rci_t* column_perms = (rci_t*) malloc(sizeof(rci_t) * n);
    CHECK_MALLOC(colmun_perms);
    CHECK_MALLOC(column_perms_copy);

    for (i = 0; i < n; i++) column_perms[i] = i;

    mzd_t* Gtemp = mzd_copy(NULL, G);

    // Ensure that we work with a systematic generator matrix
    rref_to_systematic(Gtemp, column_perms);
    mzd_t* Glw = mzd_submatrix(NULL, Gtemp, 0, n/2, n/2, n);


    // We store v=(a,b,c),k=(Z[a]+Z[b]+Z[c])
    hashtable* ht = hashtable_init(1 << sigma, 1);

    uint64_t delta1, delta2;

    for (comb1[0] = 0; comb1[0]  < n/4; comb1[0]++) {
        for (comb1[1] = 0; comb1[1] < n/4; comb1[1]++) {
            if (comb1[0] != comb1[1]) {
                delta1 = mxor(mzd_row(Glw, comb1[0]), mzd_row(Glw, comb1[1]), sigma);
                hashtable_insert(ht, &comb1, 2*sizeof(rci_t), &delta1, 8);
            }
        }
    }

    for (comb2[0] = n/4; comb2[0]  < n/2; comb2[0]++) {
        for (comb2[1] = n/4; comb2[1] < n/2; comb2[1]++) {
            if (comb2[0] != comb2[1]) {
                delta2 = mxor(mzd_row(Glw, comb2[0]), mzd_row(Glw, comb2[1]), sigma);
                elem* ret = hashtable_retrieve(ht, &delta2, 8);

                if (ret) {
                    // TODO recreate the codeword from ret->data (=comb1) and comb2
                }
            }
        }
    }


    // TODO: next iset using canteaut&chabaud

    // TODO: permute back the columns
}
