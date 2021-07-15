#ifndef COMBINATION_H
#define COMBINATION_H

typedef struct comb_t_ {
    uint16_t *combination;
    uint8_t done;
    uint64_t k;
    uint64_t n;
} comb_t;

typedef struct comb_diff_t_ {
    uint64_t to_add;
    uint64_t to_del;
} comb_diff_t;

void print_combination(comb_t *comb_struct);
void init_combination(comb_t *comb_struct, uint16_t *combination,  uint64_t k, uint64_t n);
void next_combination(comb_t *comb_struct, comb_diff_t *comb_diff);

#endif /* COMBINATION_H */
