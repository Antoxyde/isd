#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>

#include "utils.h"
#include "combinations.h"

/* Credits for this goes to https://github.com/NicsTr/binary_masking/combinations.c 
 * Which implements Knuth's revolving-door algorithm (TAOCP 7.2.1.3)
 * */

void print_combination(comb_t *comb_struct) {
    uint64_t i;
    for (i = 0; i < comb_struct->k; i++) {
        printf("%d ", comb_struct->combination[i]);
    }
    printf("\n");
}

void init_combination(comb_t *comb_struct, uint16_t *combination,  uint64_t k, uint64_t n) {
    uint64_t j;
    for (j = 0; j < k; j++) {
        combination[j] = j;
    }
    combination[k] = n;

    comb_struct->combination = combination;
    comb_struct->done = 0;
    comb_struct->k = k;
    comb_struct->n = n;
}

/*
 *  Given the current `combination` state, and the parameter `k`,
 *  modify in place `combination` to the next combination of `k` among `n`,
 *  where `n` is defined during the initialization of the internal state
 *  `combination`.
 *
 *  Return the indexes that just changed or -1 if its the last combination
 *
 */
void next_combination(comb_t *comb_struct, comb_diff_t *comb_diff) {
    uint64_t j;
    uint16_t *comb = comb_struct->combination;

    if (comb_struct->k % 2) {
        // k is odd
        // R3 : easy case
        if (comb[0] + 1 < comb[1]) {
            comb_diff->to_del = comb[0];
            comb_diff->to_add = comb[0] + 1;
            comb[0]++;
            return;
        }

        j = 1;
    } else {
        // k is even
        // R3 : easy case
        if (comb[0] > 0) {
            comb_diff->to_del = comb[0];
            comb_diff->to_add = comb[0] - 1;
            comb[0]--;
            return;
        }

        j = 1;
        // R5 : try to increase comb[j]
        // at this point comb[j-1] = j-1
        if (comb[j] + 1 < comb[j+1]) {
            comb_diff->to_del = comb[j - 1];
            comb_diff->to_add = comb[j] + 1;
            comb[j-1] = comb[j];
            comb[j]++;
            return;
        }
        j++;

    }

    while (j < comb_struct->k) {
        // R4 : try to decrease comb[j]
        // at this point comb[j] = comb[j-1] + 1
        if (comb[j] > j) {
            comb_diff->to_del = comb[j];
            comb_diff->to_add = j - 1;
            comb[j] = comb[j-1];
            comb[j-1] = j-1;
            break;
        }
        j++;

        // R5 : try to increase comb[j]
        // at this point comb[j-1] = j-1
        if (comb[j] + 1 < comb[j+1]) {
            comb_diff->to_del = comb[j - 1];
            comb_diff->to_add = comb[j] + 1;
            comb[j-1] = comb[j];
            comb[j]++;
            break;
        }
        j++;
    }
    if (j == comb_struct->k) {
        comb_struct->done = 1;
        return; // No more comb
    }
}
