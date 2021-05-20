#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>

#include "utils.h"
#include "buckets.h"

bucket** bucket_init(size_t nb_buckets) {

    bucket** buckets = calloc(nb_buckets,  sizeof(bucket*));
    if (!buckets) {
        fprintf(stderr, "Can't calloc %ld buckets, exiting\n", nb_buckets);
        exit(1);
    }

    return buckets;
}

void bucket_put(bucket** buckets,  uint64_t key, uint64_t data) {
    assert(key < (1ULL << 25));
    bucket* b = buckets[key];

    if (!b) {
        b = malloc(sizeof(bucket));
        b->maxlen = 500;
        b->curlen = 0;
        b->tab = malloc(b->maxlen * sizeof(uint64_t));
        CHECK_MALLOC(b->tab);
        buckets[key] = b;
    }
    
    if (b->curlen == b->maxlen) {
        b->maxlen += 10;
        b->tab = realloc(b->tab, b->maxlen * sizeof(uint64_t));
        CHECK_MALLOC(b->tab);
    }

    b->tab[b->curlen] = data;
    b->curlen++;
}

void bucket_free(bucket** buckets, size_t nb_buckets) {
    for (unsigned int i = 0; i < nb_buckets; i++) {
        if (buckets[i]) { 
            free(buckets[i]->tab);
            free(buckets[i]);
        }
    }

    free(buckets);
}
